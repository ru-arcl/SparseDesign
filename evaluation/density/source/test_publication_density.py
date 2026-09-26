#!/usr/bin/env python3
"""Meaningful numerical and fail-closed checks for the archived-count analyzer."""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics
import tempfile
import unittest

import publication_density as density

OUT = None


class DensityChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.products, cls.panel, cls.subset, cls.clusters, cls.provenance = density.verified_inputs(OUT)

    def test_historical_independent_audit(self):
        audit = json.loads((OUT/"verification/reference-audit.json").read_text())
        for objective in ("human-lambda0", "human-lambda4"):
            observed = density.describe(self.products[objective]); expected = audit["human"][objective]
            self.assertAlmostEqual(observed["fit"]["beta"], expected["fits_by_minimum_aa"]["0"]["beta"], places=12)
            self.assertAlmostEqual(observed["density"]["median"], expected["density_median"], places=12)
            self.assertAlmostEqual(observed["retention"]["median"], expected["retention_median"], places=12)
        for objective, expected in audit["codon_subset"].items():
            rows = [r for r in self.products[objective] if density.key(r) in self.subset]
            self.assertAlmostEqual(density.fit(rows)["beta"], expected["fits_by_minimum_aa"]["0"]["beta"], places=12)

    def test_weighted_fit_against_stdlib_regression(self):
        # Integer frequency weights must equal explicit repetition, allowing
        # comparison with an independently implemented library regression.
        rows = self.products["human-lambda0"][::97]
        weights = [(i%3)+1 for i in range(len(rows))]
        repeated = [r for r,w in zip(rows,weights) for _ in range(w)]
        expected = statistics.linear_regression([math.log(int(r["rna_nt"])) for r in repeated], [math.log(int(r["candidates"])) for r in repeated])
        observed = density.fit(rows,weights)
        self.assertAlmostEqual(observed["beta"],expected.slope,places=12)
        self.assertAlmostEqual(observed["intercept_ln"],expected.intercept,places=12)

    def test_paired_bootstrap_preserves_scale_invariance(self):
        # Multiplying all candidate counts by two changes only the intercept.
        # Every paired draw, including every cluster-weighting variant, must
        # therefore have the same slope at the two artificial endpoints.
        products = dict(self.products)
        products["human-lambda4"] = [dict(r,candidates=str(2*int(r["candidates"]))) for r in self.products["human-lambda0"]]
        draws = density.bootstrap(products,self.panel,self.clusters,8,20260716)
        for draw in draws:
            for name in ("within","cluster","equal_cluster","calibrated_equal_cluster"):
                self.assertAlmostEqual(draw[name+"_lambda0"],draw[name+"_lambda4"],places=10)

    def test_changed_csv_is_rejected(self):
        with tempfile.TemporaryDirectory(prefix="density-drift-") as directory:
            temp = Path(directory)/"density"; shutil.copytree(OUT,temp)
            csv = temp/"inputs/merged/human-lambda0.csv.gz"
            raw = bytearray(csv.read_bytes()); raw[-10] ^= 1; csv.write_bytes(raw)
            with self.assertRaisesRegex(ValueError,"hash drift"):
                density.verified_inputs(temp)

    def test_duplicate_cluster_is_rejected_beyond_hash_check(self):
        with tempfile.TemporaryDirectory(prefix="density-membership-") as directory:
            temp = Path(directory)/"density"; shutil.copytree(OUT,temp)
            path = temp/"clustering/clusters.tsv"; raw = path.read_text();path.write_text(raw+raw.splitlines()[0]+"\n")
            provenance = temp/"clustering/run-provenance.json"; metadata = json.loads(provenance.read_text())
            metadata["clusters_tsv_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            provenance.write_text(json.dumps(metadata))
            with self.assertRaisesRegex(ValueError,"Duplicate cluster membership"):
                density.verified_inputs(temp)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--out",type=Path,required=True)
    args,unittest_args=parser.parse_known_args();OUT=args.out.resolve()
    unittest.main(argv=[__file__,*unittest_args])
