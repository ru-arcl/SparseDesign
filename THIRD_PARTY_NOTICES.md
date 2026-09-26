# Third-party notices

The license for original SparseDesign work does not replace the terms or attribution of
third-party material. This package preserves the following separate notices.

## Turner 2004 numerical parameters distributed by ViennaRNA

`src/vienna/energy_parameter.h`, `intl11.h`, `intl21.h` and `intl22.h` are mechanically
generated numerical tables from ViennaRNA 2.7.2's `RNA.parameter_set_rna_turner2004`.
The same generated tables appear in frozen solver snapshots under `evaluation/`.
The generator is original project code; the upstream parameter source is separately credited.

- Upstream project: <https://www.tbi.univie.ac.at/RNA/>
- Source tag: <https://github.com/ViennaRNA/ViennaRNA/tree/v2.7.2>
- Parameter source digest: `2a43345a495850cfd2e0a78c57c6e02085e6df3c53496fe3289dc294d21732ad`
- Preserved upstream terms: [LICENSES/ViennaRNA-COPYING.txt](LICENSES/ViennaRNA-COPYING.txt),
  copied from the ViennaRNA 2.7.2 Python distribution; official tagged text is
  <https://raw.githubusercontent.com/ViennaRNA/ViennaRNA/v2.7.2/COPYING>.

Credit the ViennaRNA authors and the Institute for Theoretical Chemistry, University of
Vienna. The package reference is Lorenz et al., *ViennaRNA Package 2.0*, Algorithms for
Molecular Biology 6:26 (2011), <https://doi.org/10.1186/1748-7188-6-26>.
Thermodynamic-model reference: Mathews et al., PNAS 101:7287–7292 (2004),
<https://doi.org/10.1073/pnas.0401799101>.
ViennaRNA executables, shared libraries and its folding implementation are not bundled.

## Ward correctness fixtures

The files under `evaluation/epyc/validation/inputs/ward/` derive from
<https://github.com/maxhwardg/mrna_folding_comparison>, commit
`1391cc795309b3da82ef1ba9b21516e8327134d9`. Their original MIT license and
`Copyright (c) 2024 Max Ward` notice remain in that directory and are also copied to
[LICENSES/Ward-MIT.txt](LICENSES/Ward-MIT.txt). These files retain their MIT terms.
Reference: Ward, Richardson and Metkar, *mRNA folding algorithms for structure and codon
optimization*, Briefings in Bioinformatics 26(4):bbaf386 (2025),
<https://doi.org/10.1093/bib/bbaf386>.

## Numerical inputs and recorded results

Human/yeast codon-frequency tables and the additional evaluation tables are attributed in
`PROVENANCE.md` and their study manifests. The codon-distribution source is the
[Kazusa Codon Usage Database](https://www.kazusa.or.jp/codon/) (Nakamura, Gojobori and Ikemura,
Nucleic Acids Research 28:292, 2000; <https://doi.org/10.1093/nar/28.1.292>), with the
rounded/full-precision transformations recorded by each study.

Frozen protein panels are selected and reformatted from UniProtKB/Swiss-Prot, credited to
the UniProt Consortium. Their accessions, original release and unchanged sequence identities
are recorded in the input manifests. UniProt applies
[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/)
to copyrightable database content; see its [license notice](https://www.uniprot.org/help/license/).
The panel selection and packaging do not assign original authorship or replace these data
terms with the software's noncommercial license. Reference: The UniProt Consortium,
*UniProt: the Universal Protein Knowledgebase in 2025*, Nucleic Acids Research 53:D609–D617,
<https://doi.org/10.1093/nar/gkae1010>.

Historical comparator records contain experimental outputs, timings and source/binary
identifiers, not the comparator's implementation. The restricted LinearDesign distribution,
its binaries, shared libraries and headers are not part of this release. External baseline
software must be obtained from its own publisher under its own terms. Citing an algorithm
or retaining a numerical comparison does not assign its software license to SparseDesign.
