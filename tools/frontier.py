"""SparseDesign supported MFE/CAI-envelope and inverse-lambda tools.

This tool adaptively explores the objective tradeoff. Because the design objective
    F(r, lambda) = MFE(r) + lambda * C(r),   C(r) = sum_codons (-log w) = the CAI cost
is LINEAR in lambda, the optimal value f*(lambda) = min_r F is concave piecewise-linear, and the
weighted-sum solutions are the supported objective lines on that lower envelope.
Eisner-Severance (1976) / Aneja-Nair (1979) dichotomic search numerically discovers the extreme
envelope segments with about two exact solves per representative instead of a blind grid. It does
not claim to enumerate unsupported nondominated designs, and one representative is retained when
multiple sequences define the same objective line.

The inverse tool takes a protein and one synonymous mRNA, folds that fixed sequence with the same
Turner model, and numerically infers the interval of nonnegative lambda values for which its
objective line lies on the discovered supported lower envelope. The answer can be an interval, one
breakpoint, or empty; lambda is not generally identifiable as one scalar from a sequence alone.

Two numerical policies govern this calculation:
  * Derive C from the designed sequence's full-precision codon weights, not the binary's printed
    3-decimal CAI (that coarse C corrupts the intersection-lambda geometry and misses shallow vertices).
  * the left endpoint uses lambda=0+ (a tiny positive lambda) as a numerical policy that favors
    higher CAI among near-tied min-MFE results; it is not a lexicographic exactness guarantee
    -- avoiding the dominated point a literal lambda=0 solve returns.
"""
import hashlib
import math
import os
import shutil
import subprocess
import collections

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_BIN = os.environ.get("SPARSEDESIGN_BIN",
                      os.environ.get("LINEARDESIGN_CLEAN_BIN",
                                     os.path.join(_REPO, "sparsedesign")))
_TABLE = os.environ.get("CODON_TABLE",
                        os.path.join(_REPO, "data",
                                     "codon_usage_freq_table_human.csv"))
_THREADS = os.environ.get("DESIGN_THREADS", "1")
_TIMEOUT = int(os.environ.get("DESIGN_TIMEOUT", "1800"))

_LAM_MIN = 1e-4     # numerical lambda=0+ endpoint; not an exact lexicographic tie-break
_LAM_MAX = 50.0     # initial upper search bound; not assumed to be the max-CAI endpoint
_EPS = 1e-7         # "strictly below the chord" tolerance (float noise floor; MFE granularity is 0.01)
_C_EPS = 1e-12
_MAX_RECURSION_DEPTH = 100
_HARD_LAM_MAX = 1e12
_TIE_BREAK_ENERGY_MARGIN = 0.004
_INTERVAL_RELATIVE_EPS = 1e-10

_FRONTIER_KIND = "epsilon-supported-envelope-segment-representatives"
_OBJECTIVE = "minimize MFE_kcal_mol + lambda * sum_codons(-log(relative_adaptiveness))"
_SEMANTICS = (
    "Numerically discovered extreme weighted-sum lower-envelope segments. Unsupported "
    "nondominated designs and distinct designs supported only at a zero-width tie breakpoint are "
    "not enumerated; tied sequences with the same MFE/codon-cost line share one representative."
)


def _sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve_executable(path):
    if os.path.dirname(path):
        return os.path.realpath(os.path.abspath(path))
    resolved = shutil.which(path)
    if resolved:
        return os.path.realpath(resolved)
    return os.path.abspath(path)


def backend_identity(bin_path=_BIN, table=_TABLE, threads=_THREADS,
                     timeout=_TIMEOUT, effective_threads=None):
    """Return auditable identity for the executable and codon table actually configured."""
    binary_path = _resolve_executable(bin_path)
    table_path = os.path.realpath(os.path.abspath(table))
    requested_threads = int(threads)
    return {
        "id": os.path.basename(binary_path),
        "binary_path": binary_path,
        "binary_sha256": _sha256_file(binary_path),
        "codon_table_path": table_path,
        "codon_table_sha256": _sha256_file(table_path),
        "threads": requested_threads,
        "requested_threads": requested_threads,
        "effective_threads": effective_threads,
        "timeout_seconds": timeout,
    }


def _single_fasta_sequence(raw, label):
    """Normalize one raw/FASTA sequence and reject accidental record concatenation."""
    nonempty = [line.strip() for line in (raw or "").splitlines() if line.strip()]
    headers = [index for index, line in enumerate(nonempty) if line.startswith(">")]
    if len(headers) > 1:
        raise ValueError(
            "multiple FASTA %s records are not supported by this single-job API" % label
        )
    if headers and headers[0] != 0:
        raise ValueError("a FASTA %s header must precede its sequence" % label)
    sequence_lines = nonempty[1:] if headers else nonempty
    sequence = "".join("".join(line.split()) for line in sequence_lines).upper()
    if not sequence:
        raise ValueError("%s input is empty" % label)
    return sequence


def _load_weight_model(table_path):
    """Parse codon weights with the same amino-acid and frequency semantics as C++.

    Returns ``(codon_weights, minimum_cost_by_amino_acid)``. The C++ parser treats the first
    character of a nonempty amino-acid field as its one-letter code, so this adapter must do the
    same: otherwise a custom table mixing fields such as ``A`` and ``Ala`` would be optimized under
    one set of relative-adaptiveness weights and analyzed here under another.
    """
    rows, maxf = [], collections.defaultdict(float)
    with open(table_path) as fh:
        for line in fh:
            line = line.strip().lstrip("﻿")
            if not line or line.startswith("#"):
                continue
            parts = [x.strip() for x in line.split(",")]
            if len(parts) < 3 or len(parts[0]) != 3 or not parts[1]:
                continue
            try:
                freq = float(parts[2])
            except ValueError:
                raise ValueError(
                    "invalid codon frequency for %s: %r" % (parts[0], parts[2]))
            if not math.isfinite(freq) or freq < 0:
                raise ValueError("invalid codon frequency for %s: %r" % (parts[0], parts[2]))
            codon = parts[0].upper().replace("T", "U")
            aa = parts[1][0]
            rows.append((codon, aa, freq))
            maxf[aa] = max(maxf[aa], freq)
    if not rows:
        raise ValueError("no codon rows parsed from: %s" % table_path)

    weights = {}
    minimum_cost_by_aa = {}
    for codon, aa, freq in rows:
        weight = freq / maxf[aa] if maxf[aa] > 0 else 0.0
        weight = weight if weight > 0 else 1e-9
        weights[codon] = weight
        cost = -math.log(weight)
        if aa not in minimum_cost_by_aa or cost < minimum_cost_by_aa[aa]:
            minimum_cost_by_aa[aa] = cost
    return weights, minimum_cost_by_aa


def _load_weights(table_path):
    """Compatibility helper returning codon -> relative adaptiveness."""
    return _load_weight_model(table_path)[0]


class _Solver:
    def __init__(self, protein, bin_path=_BIN, table=_TABLE, threads=_THREADS,
                 timeout=_TIMEOUT, on_solve=None):
        self.protein_seq = _single_fasta_sequence(protein, "protein")
        self.protein = ">p\n" + self.protein_seq + "\n"
        self.requested_threads = int(threads)
        self.threads = str(self.requested_threads)
        self.timeout, self.on_solve = timeout, on_solve
        self.effective_threads = None
        self._backend = backend_identity(
            bin_path, table, threads=self.requested_threads, timeout=timeout)
        self.bin = self._backend["binary_path"]
        self.table = self._backend["codon_table_path"]
        self.w, minimum_cost_by_aa = _load_weight_model(self.table)
        missing = sorted(set(self.protein_seq) - set(minimum_cost_by_aa))
        if missing:
            raise ValueError(
                "codon-usage table has no codons for amino acid(s): %s"
                % ", ".join(missing))
        self.minimum_codon_cost = math.fsum(
            minimum_cost_by_aa[aa] for aa in self.protein_seq)
        self.count = 0

    def _cai_cost(self, seq):
        if len(seq) % 3:
            raise ValueError("mRNA length is not divisible by three")
        costs = []
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon not in self.w or self.w[codon] <= 0:
                raise ValueError("codon is absent from the codon-usage table: %s" % codon)
            costs.append(-math.log(self.w[codon]))
        return math.fsum(costs)

    def _point_from_output(self, out, lam):
        mfe = float(next(line for line in out.splitlines() if "free energy" in line)
                    .split("free energy:")[1].split("kcal")[0])
        seq = next(line for line in out.splitlines() if "mRNA sequence" in line).split(":", 1)[1].strip()
        struct = next(line for line in out.splitlines() if "mRNA structure" in line).split(":", 1)[1].strip()
        thread_lines = [line for line in out.splitlines() if line.startswith("Threads:")]
        if thread_lines:
            self.effective_threads = int(thread_lines[0].split(":", 1)[1].split()[0])
        cost = self._cai_cost(seq)
        cai = math.exp(-cost / (len(seq) // 3)) if seq else 1.0
        return {"lambda": lam, "mfe": mfe, "cai": cai, "C": cost,
                "sequence": seq, "structure": struct}

    def _run(self, args):
        try:
            proc = subprocess.run(args, input=self.protein, capture_output=True, text=True,
                                  timeout=self.timeout)
        except subprocess.TimeoutExpired:
            raise RuntimeError("optimizer exceeded the %g-second per-solve timeout" % self.timeout)
        if proc.returncode:
            detail = proc.stderr.strip() or proc.stdout.strip() or "optimizer failed"
            raise RuntimeError(detail)
        return proc.stdout

    def solve(self, lam):
        """Run the exact optimizer at tradeoff lambda; return one frontier point."""
        self.count += 1
        out = self._run([
            self.bin, "-l", str(lam), "-j", self.threads, "-c", self.table, "-v"])
        point = self._point_from_output(out, lam)
        if self.on_solve:
            self.on_solve(self.count, "optimize", lam)
        return point

    def evaluate(self, mrna):
        """Fold and score one supplied synonymous mRNA through the production binary."""
        seq = _single_fasta_sequence(mrna, "mRNA").replace("T", "U")
        self.count += 1
        out = self._run([self.bin, "--evaluate-mrna", seq, "-j", self.threads,
                         "-c", self.table, "-v"])
        point = self._point_from_output(out, None)
        if self.on_solve:
            self.on_solve(self.count, "evaluate", None)
        return point

    def tie_break_lambda(self):
        """A safe 0+ value: preserve integer-cent MFE, then prefer lower codon cost."""
        codons = len(self.protein_seq)
        max_per_codon = max((-math.log(w) for w in self.w.values() if w > 0), default=0.0)
        max_total_cost = codons * max_per_codon
        return (min(_LAM_MIN, _TIE_BREAK_ENERGY_MARGIN / max_total_cost)
                if max_total_cost else _LAM_MIN)

    def minimum_attainable_codon_cost(self):
        """Return the exact minimum C attainable for this protein under the loaded table."""
        return self.minimum_codon_cost

    def backend_identity(self):
        current = backend_identity(
            self.bin, self.table, threads=self.requested_threads,
            timeout=self.timeout, effective_threads=self.effective_threads)
        for key in ("binary_sha256", "codon_table_sha256"):
            if current[key] != self._backend[key]:
                raise RuntimeError(
                    "optimizer provenance changed during the frontier job (%s)" % key)
        return current


def _frontier_with_solver(s, lam_min, lam_max, lo=None, hi=None, objective_epsilon=_EPS,
                          codon_cost_epsilon=_C_EPS,
                          max_recursion_depth=_MAX_RECURSION_DEPTH):
    if not (0.0 < lam_min < lam_max and math.isfinite(lam_max)):
        raise ValueError("frontier bounds must satisfy 0 < lambda_min < lambda_max < infinity")
    if not (math.isfinite(objective_epsilon) and objective_epsilon >= 0.0):
        raise ValueError("objective_epsilon must be finite and nonnegative")
    if not (math.isfinite(codon_cost_epsilon) and codon_cost_epsilon >= 0.0):
        raise ValueError("codon_cost_epsilon must be finite and nonnegative")
    if not isinstance(max_recursion_depth, int) or max_recursion_depth < 0:
        raise ValueError("max_recursion_depth must be a nonnegative integer")
    lo = lo or s.solve(lam_min)
    hi = hi or s.solve(lam_max)
    verts = {lo["sequence"]: lo, hi["sequence"]: hi}

    def recurse(a, b, depth=0):
        # ``<=`` is required when callers request exact zero tolerance: identical slopes must still
        # terminate before the intersection formula divides by zero.
        if abs(a["C"] - b["C"]) <= codon_cost_epsilon:
            return
        lam = (b["mfe"] - a["mfe"]) / (a["C"] - b["C"])
        if not (lam_min < lam < lam_max):
            return
        if depth >= max_recursion_depth:
            raise RuntimeError(
                "supported-frontier recursion depth limit (%d) reached before convergence"
                % max_recursion_depth
            )
        d = s.solve(lam)
        if (d["mfe"] + lam * d["C"]
                < a["mfe"] + lam * a["C"] - objective_epsilon
                and d["sequence"] not in (a["sequence"], b["sequence"])):
            verts[d["sequence"]] = d
            recurse(a, d, depth + 1)
            recurse(d, b, depth + 1)

    recurse(lo, hi)
    pts = list(verts.values())
    front = [x for x in pts if not any(
        y["mfe"] <= x["mfe"] + objective_epsilon
        and y["C"] <= x["C"] + codon_cost_epsilon
        and (y["mfe"] < x["mfe"] - objective_epsilon
             or y["C"] < x["C"] - codon_cost_epsilon)
        for y in pts)]
    # math.fsum makes identical codon-cost multisets stable; keep one representative for an exactly
    # identical floating objective line without merging merely-nearby slopes at large lambda.
    unique = {}
    for point in front:
        unique.setdefault((point["mfe"], point["C"]), point)
    return sorted(unique.values(), key=lambda x: (x["mfe"], -x["C"]))


def supported_frontier(protein, lambda_min=None, lambda_max=_LAM_MAX,
                       hard_lambda_max=_HARD_LAM_MAX, objective_epsilon=_EPS,
                       codon_cost_epsilon=_C_EPS,
                       max_recursion_depth=_MAX_RECURSION_DEPTH, **kw):
    """Return a provenance-bearing numerical supported MFE/CAI envelope.

    ``lambda_max`` is doubled as needed until the exact minimum attainable codon-cost endpoint is
    observed. That minimum is normally zero, but can be positive for a custom table whose entire
    synonymous group has zero frequencies and is therefore floored by the optimizer.
    ``lambda_min=None`` selects a protein-length-safe numerical 0+ tie-break value. Search
    convergence is subject to the reported chord and codon-cost tolerances; ``complete`` is
    therefore deliberately false rather than a mathematical exactness claim.
    """
    s = _Solver(protein, **kw)
    if lambda_min is not None and not (math.isfinite(lambda_min) and lambda_min > 0.0):
        raise ValueError("lambda_min must be finite and positive, or None for automatic 0+")
    if not (math.isfinite(lambda_max) and lambda_max > 0.0):
        raise ValueError("lambda_max must be finite and positive")
    if not (math.isfinite(hard_lambda_max) and hard_lambda_max >= lambda_max):
        raise ValueError("hard_lambda_max must be finite and at least lambda_max")

    effective_min = s.tie_break_lambda() if lambda_min is None else lambda_min
    if not effective_min < lambda_max:
        raise ValueError("effective lambda_min must be less than lambda_max")
    lo = s.solve(effective_min)
    minimum_codon_cost = s.minimum_attainable_codon_cost()
    effective_max = lambda_max
    hi = s.solve(effective_max)
    while hi["C"] != minimum_codon_cost:
        if effective_max >= hard_lambda_max:
            raise RuntimeError(
                "failed to reach the minimum-codon-cost endpoint (C=%.12g) by lambda=%g"
                % (minimum_codon_cost, hard_lambda_max)
            )
        effective_max = min(hard_lambda_max, effective_max * 2.0)
        hi = s.solve(effective_max)

    effective_codon_epsilon = min(
        codon_cost_epsilon,
        objective_epsilon / effective_max if objective_epsilon > 0.0 else 0.0,
    )
    front = _frontier_with_solver(
        s, effective_min, effective_max, lo=lo, hi=hi,
        objective_epsilon=objective_epsilon,
        codon_cost_epsilon=effective_codon_epsilon,
        max_recursion_depth=max_recursion_depth,
    )
    zero_cost_endpoint_reached = hi["C"] == 0.0
    minimum_cost_endpoint_reached = hi["C"] == minimum_codon_cost
    return {
        "schema": "sparsedesign.supported-frontier",
        "schema_version": 1,
        "frontier": front,
        "frontier_kind": _FRONTIER_KIND,
        "objective": _OBJECTIVE,
        "semantics": _SEMANTICS,
        "configured_lambda_bounds": {
            "lambda_min": lambda_min,
            "lambda_max": lambda_max,
            "hard_lambda_max": hard_lambda_max,
            "lambda_units": "kcal/mol per unit -ln(relative_adaptiveness)",
        },
        "effective_lambda_bounds": {
            "lambda_min": effective_min,
            "lambda_max": effective_max,
        },
        "zero_plus_lambda": effective_min if lambda_min is None else None,
        "numerical_lower_bound": effective_min,
        "search_policy": {
            "lower_endpoint": ("length-safe-zero-plus" if lambda_min is None
                               else "configured-positive-bound"),
            "upper_endpoint": "double-until-minimum-attainable-codon-cost",
        },
        "tolerances": {
            "objective_epsilon": objective_epsilon,
            "configured_codon_cost_epsilon": codon_cost_epsilon,
            "effective_codon_cost_epsilon": effective_codon_epsilon,
            "max_recursion_depth": max_recursion_depth,
            "zero_plus_energy_margin_kcal_mol": _TIE_BREAK_ENERGY_MARGIN,
        },
        "minimum_attainable_codon_cost": minimum_codon_cost,
        "minimum_cost_endpoint_reached": minimum_cost_endpoint_reached,
        "zero_cost_endpoint_reached": zero_cost_endpoint_reached,
        "upper_endpoint_codon_cost": hi["C"],
        "solve_count": s.count,
        "complete": False,
        "complete_within_effective_bounds": False,
        "search_converged": True,
        "epsilon_complete_within_effective_bounds": True,
        "completeness_scope":
            "epsilon-filtered-extreme-supported-envelope-segments-within-effective-bounds",
        "point_lambda_semantics":
            "witness solve value, not a support interval or guaranteed breakpoint",
        "backend": s.backend_identity(),
    }


def frontier(protein, lambda_min=_LAM_MIN, lambda_max=_LAM_MAX, **kw):
    """Compatibility wrapper returning supported designs and solve count as a tuple.

    The legacy explicitly bounded behavior is retained. New callers should use
    :func:`supported_frontier` for length-safe 0+, endpoint completion, bounds, tolerances, and
    backend provenance.
    """
    s = _Solver(protein, **kw)
    front = _frontier_with_solver(s, lambda_min, lambda_max)
    return front, s.count


def infer_lambda(protein, mrna, initial_lambda_max=_LAM_MAX,
                 hard_lambda_max=_HARD_LAM_MAX, **kw):
    """Numerically infer a nonnegative lambda interval supporting ``mrna`` for ``protein``.

    The upper search endpoint is doubled until the optimizer reaches the minimum codon cost
    attainable for this protein and codon table; that cost is normally zero but need not be for
    custom all-zero synonymous groups. The lower endpoint is a length-safe 0+. The conclusion is
    subject to the reported numerical tolerances; ``supported`` is false when the target line does
    not touch the discovered supported envelope.
    """
    if not (math.isfinite(initial_lambda_max) and initial_lambda_max > 0):
        raise ValueError("initial_lambda_max must be finite and positive")
    if not (math.isfinite(hard_lambda_max) and hard_lambda_max >= initial_lambda_max):
        raise ValueError("hard_lambda_max must be finite and at least initial_lambda_max")

    s = _Solver(protein, **kw)
    target = s.evaluate(mrna)
    minimum_codon_cost = s.minimum_attainable_codon_cost()
    lam_min = s.tie_break_lambda()
    lo = s.solve(lam_min)
    lam_max = initial_lambda_max
    hi = s.solve(lam_max)
    while hi["C"] != minimum_codon_cost:
        if lam_max >= hard_lambda_max:
            raise RuntimeError(
                "failed to reach the minimum-codon-cost endpoint (C=%.12g) by lambda=%g"
                % (minimum_codon_cost, hard_lambda_max))
        lam_max = min(hard_lambda_max, lam_max * 2.0)
        hi = s.solve(lam_max)

    effective_codon_epsilon = min(
        _C_EPS, _EPS / lam_max if _EPS > 0.0 else 0.0)
    front = _frontier_with_solver(
        s, lam_min, lam_max, lo=lo, hi=hi,
        codon_cost_epsilon=effective_codon_epsilon)
    lower, upper = 0.0, math.inf
    impossible = False
    for competitor in front:
        delta_c = target["C"] - competitor["C"]
        delta_e = competitor["mfe"] - target["mfe"]
        if abs(delta_c) <= effective_codon_epsilon:
            if target["mfe"] > competitor["mfe"] + _EPS:
                impossible = True
                break
        elif delta_c > 0:
            upper = min(upper, delta_e / delta_c)
        else:
            lower = max(lower, delta_e / delta_c)

    lower = max(0.0, lower)
    scale = max(1.0, abs(lower), abs(upper) if math.isfinite(upper) else 1.0)
    supported = (not impossible and upper >= -_EPS
                 and lower <= upper + _INTERVAL_RELATIVE_EPS * scale)
    if supported:
        upper = max(0.0, upper) if math.isfinite(upper) else upper
        if (math.isfinite(upper)
                and abs(upper - lower) <= _INTERVAL_RELATIVE_EPS * scale):
            lower = upper = (lower + upper) / 2.0
        representative = ((lower + upper) / 2.0 if math.isfinite(upper)
                          else max(lower, lam_max))
        unbounded = math.isinf(upper)
        reported_upper = None if unbounded else upper
    else:
        lower = upper = representative = None
        unbounded = False
        reported_upper = None

    return {"schema": "sparsedesign.lambda-support", "schema_version": 1,
            "supported": supported, "lambda_min": lower,
            "lambda_max": reported_upper, "unbounded": unbounded,
            "representative_lambda": representative, "target": target,
            "frontier_vertices": len(front), "solve_count": s.count,
            "search_lambda_max": lam_max, "frontier_kind": _FRONTIER_KIND,
            "objective": _OBJECTIVE, "semantics": _SEMANTICS,
            "configured_lambda_bounds": {
                "lambda_min": 0.0,
                "initial_lambda_max": initial_lambda_max,
                "hard_lambda_max": hard_lambda_max,
                "lambda_units": "kcal/mol per unit -ln(relative_adaptiveness)",
            },
            "effective_lambda_bounds": {
                "lambda_min": lam_min,
                "lambda_max": lam_max,
            },
            "zero_plus_lambda": lam_min,
            "search_policy": {
                "lower_endpoint": "length-safe-zero-plus",
                "upper_endpoint": "double-until-minimum-attainable-codon-cost",
            },
            "tolerances": {
                "objective_epsilon": _EPS,
                "configured_codon_cost_epsilon": _C_EPS,
                "effective_codon_cost_epsilon": effective_codon_epsilon,
                "max_recursion_depth": _MAX_RECURSION_DEPTH,
                "support_interval_relative_epsilon": _INTERVAL_RELATIVE_EPS,
                "zero_plus_energy_margin_kcal_mol": _TIE_BREAK_ENERGY_MARGIN,
            },
            "minimum_attainable_codon_cost": minimum_codon_cost,
            "minimum_cost_endpoint_reached": hi["C"] == minimum_codon_cost,
            "zero_cost_endpoint_reached": hi["C"] == 0.0,
            "upper_endpoint_codon_cost": hi["C"],
            "complete": False,
            "complete_within_effective_bounds": False,
            "search_converged": True,
            "epsilon_complete_within_effective_bounds": True,
            "completeness_scope":
                "epsilon-filtered-lambda-support-over-discovered-extreme-envelope-segments",
            "point_lambda_semantics":
                "witness solve value, not a support interval or guaranteed breakpoint",
            "backend": s.backend_identity()}


def _read_fasta_sequence(path):
    with open(path) as handle:
        return _single_fasta_sequence(handle.read(), "mRNA")


if __name__ == "__main__":
    import argparse
    import json
    import sys
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--infer-mrna", metavar="FASTA",
                        help="report the lambda interval supporting this mRNA for stdin protein")
    parser.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    parser.add_argument("-j", "--threads", default=_THREADS)
    parser.add_argument("-c", "--codon-table", default=_TABLE)
    parser.add_argument("--binary", default=_BIN)
    args = parser.parse_args()
    prot = sys.stdin.read()
    options = {"bin_path": args.binary, "table": args.codon_table, "threads": args.threads}
    if args.infer_mrna:
        result = infer_lambda(prot, _read_fasta_sequence(args.infer_mrna), **options)
        if args.json:
            print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
        else:
            backend = result["backend"]
            print("# backend %s; binary sha256 %s; codon-table sha256 %s" %
                  (backend["id"], backend["binary_sha256"],
                   backend["codon_table_sha256"]))
            target = result["target"]
            print("# target MFE %.2f kcal/mol; CAI %.6f; C %.12g" %
                  (target["mfe"], target["cai"], target["C"]))
            print("# %d discovered envelope-segment representatives; %d exact solves" %
                  (result["frontier_vertices"], result["solve_count"]))
            if result["supported"]:
                hi = ("infinity" if result["unbounded"]
                      else "%.12g" % result["lambda_max"])
                print("supporting lambda interval: [%.12g, %s]" %
                      (result["lambda_min"], hi))
                print("representative lambda: %.12g" % result["representative_lambda"])
            else:
                print("supporting lambda interval: none")
    else:
        result = supported_frontier(prot, **options)
        front, n = result["frontier"], result["solve_count"]
        if args.json:
            print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
        else:
            backend = result["backend"]
            print("# backend %s; binary sha256 %s; codon-table sha256 %s" %
                  (backend["id"], backend["binary_sha256"],
                   backend["codon_table_sha256"]))
            print("# numerical supported MFE/CAI envelope: %d segment representatives "
                  "in %d solves" % (len(front), n))
            bounds = result["effective_lambda_bounds"]
            print("# effective lambda interval [%g, %g]; epsilon-converged=%s; "
                  "minimum-cost endpoint=%s (C=%.12g)" %
                  (bounds["lambda_min"], bounds["lambda_max"],
                   str(result["search_converged"]).lower(),
                   str(result["minimum_cost_endpoint_reached"]).lower(),
                   result["minimum_attainable_codon_cost"]))
            print(f"# {'MFE':>9} {'CAI':>7}  lambda*")
            for d in front:
                print(f"  {d['mfe']:9.2f} {d['cai']:7.4f}  {d['lambda']:.4g}")
