# Source and data provenance

Original SparseDesign software and documentation are supplied under the
[Rutgers Non-commercial Research License (RU-NCRL)](LICENSE.md).
The text reproduces the user-supplied [Word document](LICENSE.docx) word for word;
only paragraph spacing and bullet presentation are adapted for Markdown.
[NOTICE](NOTICE) reproduces the Rutgers copyright notice specified in that license.
Separate terms and attribution for third-party material are listed in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Implementation lineage

The solver, numerical-table generator and deterministic test suite come from the project's
`lineardesign-clean/` development tree, with publication-audit fixes. SparseDesign is the release
name; `sparsedesign` is the canonical executable, and the build preserves `lineardesign-clean`
as a compatibility name. Internal `ldclean` identifiers and historical artifact names remain
traceable to that development tree. The project records this implementation
as independently written from published algorithmic and energy-model descriptions.
`tools/frontier.py` comes from the project's own `webservice/frontier.py`, adapted to a standalone
solver root; it uses the Python standard library. The generated parameter headers have the
specific upstream provenance below.

No restricted LinearDesign executable, shared library, header or implementation source is
included. The optional historical test driver requiring its shared library is excluded.
Frozen source snapshots accompanying evaluations identify the version actually tested;
they do not silently inherit later CLI, executable-name, or documentation changes. Renaming
the release does not relabel recorded binaries or establish new timing measurements.

## Energy parameters

The four `src/vienna/*.h` files are generated from ViennaRNA 2.7.2
`RNA.parameter_set_rna_turner2004`, whose complete source-text SHA-256 is:

```
2a43345a495850cfd2e0a78c57c6e02085e6df3c53496fe3289dc294d21732ad
```

`tools/generate_vienna_parameters.py` documents index conversion and unknown-base filling,
enforces this digest, and compares every generated header. Upstream ViennaRNA terms are
preserved in `LICENSES/ViennaRNA-COPYING.txt`. ViennaRNA itself is neither bundled nor needed
for a normal solver build; it is an optional independent validation dependency.

## Input and result identity

The bundled human and yeast codon tables are the numerical inputs used by the project.
Their source-snapshot hashes are in `SOURCE_SHA256SUMS`; campaign manifests identify each
actual table used. The density study preserves the original host-table metadata, source
counts, rounded/full-precision conventions, UniProt sequence identities and input hashes.
Ward regression fixtures retain their original MIT notice.

Completed evaluations have their own frozen inputs, result summaries, source identities and
checksums under `evaluation/`. Historical measurements use older binaries and the hardware
and protocol named by each dataset; they are not timings of a fresh build of this release.
`evaluation/historical/` contains compact records and a relocatable checker, with no external
solver implementation or executable. The [benchmark reference](BENCHMARKS.md) summarizes
the current [primary EPYC studies](evaluation/epyc/README.md) and the separate i9 Dp427c
commodity-PC demonstration in `evaluation/commodity-pc/`. Other superseded i9 campaigns
are omitted from this release; their original records remain in the development repository.
The active solver and codon tables retain the measured numerical
implementation; differences from the EPYC source snapshots are existing branding comments
and the CLI help banner. The benchmark update requires no solver changes.

## Releasing and citing a snapshot

`CITATION.cff` records the SparseDesign name, authors, license and
[source repository](https://github.com/ru-arcl/SparseDesign). Cite the actual commit or release
identifier together with the manuscript; no DOI is assigned here.
`SOURCE_SHA256SUMS` covers the top-level release
instructions/notices and the active solver, test, data and tool files. Evaluation studies use
separate manifests so a historical result keeps the identity of its original source snapshot.
Check the active tree with `python3 tools/update_source_manifest.py --check`; regenerate this
manifest only after reviewing intentional source or documentation changes. Preserve measured-record
and source-snapshot identities; refresh a containing evaluation manifest when reviewed documentation
changes alter its covered files.
