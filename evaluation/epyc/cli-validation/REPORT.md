# Independent checks of controlled CLI outputs

Independent translation and returned-structure energy/admissibility checks only. No fixed-sequence refolding, new sequence design, ensemble claim or independent optimality proof.

Original study complete: False. CLI outcomes: 70/70; statuses: {'ok': 70}.
Checked 70 successful executions and 5 distinct RNA/structure pairs; every successful output passes: True.

The frozen raw-record verifier runs first. Translation uses the standard genetic code; codon penalty and CAI are recomputed from the original hashed table. Structures must have canonical pairs, hairpins of at least three unpaired bases, and at most 30 unpaired bases in each bulge/internal loop.

ViennaRNA 2.7.2 loads Turner 2004 at 37 C, dangles 0, special hairpins enabled and salt 1.021 M. OPTION_EVAL_ONLY must leave both DP matrix pointers empty. eval_structure_pt returns integer centikcal/mol. Printed energies must match exactly; lambda-0 DP energies must equal that integer and realized diagnostics must equal the exact source divide/multiply round trip. Weighted objectives allow 1e-4 centikcal/mol numerical tolerance; codon penalty and CAI tolerances are 1e-8 and 1e-12. All individual checks and differences remain in JSON.

| Phase | Accession | Lambda | Threads | Layout | Repeat | Status | Checks pass |
|---|---|---:|---:|---|---:|---|---|
| scaling | Q54TT4 | 4 | 1 | square | 2 | ok | True |
| scaling | Q54TT4 | 0 | 4 | square | 1 | ok | True |
| scaling | Q54TT4 | 0 | 4 | square | 3 | ok | True |
| scaling | Q54TT4 | 0 | 16 | square | 1 | ok | True |
| scaling | Q54TT4 | 0 | 2 | square | 1 | ok | True |
| scaling | Q54TT4 | 4 | 2 | square | 1 | ok | True |
| scaling | Q54TT4 | 0 | 2 | square | 2 | ok | True |
| scaling | Q54TT4 | 4 | 2 | square | 3 | ok | True |
| scaling | Q54TT4 | 0 | 8 | square | 1 | ok | True |
| scaling | Q54TT4 | 4 | 16 | square | 2 | ok | True |
| scaling | Q54TT4 | 4 | 4 | square | 3 | ok | True |
| scaling | Q54TT4 | 4 | 1 | square | 3 | ok | True |
| scaling | Q54TT4 | 4 | 8 | square | 2 | ok | True |
| scaling | Q54TT4 | 4 | 16 | square | 1 | ok | True |
| scaling | Q54TT4 | 0 | 1 | square | 3 | ok | True |
| scaling | Q54TT4 | 0 | 16 | square | 2 | ok | True |
| scaling | Q54TT4 | 4 | 16 | square | 3 | ok | True |
| scaling | Q54TT4 | 0 | 4 | square | 2 | ok | True |
| scaling | Q54TT4 | 4 | 4 | square | 2 | ok | True |
| scaling | Q54TT4 | 4 | 1 | square | 1 | ok | True |
| scaling | Q54TT4 | 4 | 8 | square | 3 | ok | True |
| scaling | Q54TT4 | 4 | 4 | square | 1 | ok | True |
| scaling | Q54TT4 | 0 | 1 | square | 1 | ok | True |
| scaling | Q54TT4 | 0 | 2 | square | 3 | ok | True |
| scaling | Q54TT4 | 0 | 8 | square | 3 | ok | True |
| scaling | Q54TT4 | 0 | 8 | square | 2 | ok | True |
| scaling | Q54TT4 | 0 | 1 | square | 2 | ok | True |
| scaling | Q54TT4 | 0 | 16 | square | 3 | ok | True |
| scaling | Q54TT4 | 4 | 8 | square | 1 | ok | True |
| scaling | Q54TT4 | 4 | 2 | square | 2 | ok | True |
| scaling | Q61879 | 0 | 16 | square | 1 | ok | True |
| scaling | Q61879 | 0 | 2 | square | 1 | ok | True |
| scaling | Q61879 | 4 | 4 | square | 2 | ok | True |
| scaling | Q61879 | 4 | 1 | square | 2 | ok | True |
| scaling | Q61879 | 0 | 8 | square | 1 | ok | True |
| scaling | Q61879 | 0 | 4 | square | 2 | ok | True |
| scaling | Q61879 | 0 | 8 | square | 3 | ok | True |
| scaling | Q61879 | 4 | 16 | square | 3 | ok | True |
| scaling | Q61879 | 4 | 8 | square | 2 | ok | True |
| scaling | Q61879 | 0 | 2 | square | 2 | ok | True |
| scaling | Q61879 | 0 | 8 | square | 2 | ok | True |
| scaling | Q61879 | 0 | 4 | square | 1 | ok | True |
| scaling | Q61879 | 4 | 4 | square | 1 | ok | True |
| scaling | Q61879 | 0 | 16 | square | 3 | ok | True |
| scaling | Q61879 | 4 | 2 | square | 2 | ok | True |
| scaling | Q61879 | 0 | 1 | square | 1 | ok | True |
| scaling | Q61879 | 0 | 1 | square | 3 | ok | True |
| scaling | Q61879 | 4 | 1 | square | 1 | ok | True |
| scaling | Q61879 | 0 | 16 | square | 2 | ok | True |
| scaling | Q61879 | 4 | 2 | square | 1 | ok | True |
| scaling | Q61879 | 4 | 16 | square | 2 | ok | True |
| scaling | Q61879 | 4 | 8 | square | 3 | ok | True |
| scaling | Q61879 | 4 | 1 | square | 3 | ok | True |
| scaling | Q61879 | 4 | 8 | square | 1 | ok | True |
| scaling | Q61879 | 0 | 1 | square | 2 | ok | True |
| scaling | Q61879 | 0 | 2 | square | 3 | ok | True |
| scaling | Q61879 | 4 | 4 | square | 3 | ok | True |
| scaling | Q61879 | 4 | 16 | square | 1 | ok | True |
| scaling | Q61879 | 4 | 2 | square | 3 | ok | True |
| scaling | Q61879 | 0 | 4 | square | 3 | ok | True |
| workstation | NP_000100.3 | 0 | 16 | square | 1 | ok | True |
| workstation | NP_000100.3 | 0 | 16 | packed | 1 | ok | True |
| workstation | NP_000100.3 | 0 | 16 | packed | 2 | ok | True |
| workstation | NP_000100.3 | 0 | 16 | square | 2 | ok | True |
| workstation | NP_000100.3 | 0 | 16 | square | 3 | ok | True |
| workstation | NP_000100.3 | 0 | 16 | packed | 3 | ok | True |
| workstation | NP_000100.3 | 0 | 16 | packed | 4 | ok | True |
| workstation | NP_000100.3 | 0 | 16 | square | 4 | ok | True |
| workstation | NP_000100.3 | 0 | 16 | square | 5 | ok | True |
| workstation | NP_000100.3 | 0 | 16 | packed | 5 | ok | True |
