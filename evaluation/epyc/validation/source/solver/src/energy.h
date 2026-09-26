// energy.h — Turner 2004 nearest-neighbor RNA free-energy model, dangle model 0.
// Clean-room: formulas and parameters are the open-source ViennaRNA / LinearFold energy
// model (parameter tables in src/vienna/*.h). The LinearDesign paper specifies exactly this:
// "thermodynamic parameters follow LinearFold and Vienna RNAfold, except dangling ends,
// which do not contribute stability" (= dangle model 0). Energies are integers in units of
// 0.01 kcal/mol. Nucleotide encoding (Vienna): A=1, C=2, G=3, U=4.
#ifndef LDCLEAN_ENERGY_H
#define LDCLEAN_ENERGY_H

#include <string>
#include <utility>
#include <vector>

namespace ldclean {
namespace energy {

void init();                       // compatibility no-op; energy tables are immutable
bool pairable(int a, int b);       // is (a,b) a valid base pair?

// ViennaRNA pair-type code for nucleotides (a,b): 0:none 1:CG 2:GC 3:GU 4:UG 5:AU 6:UA.
int pair_type(int a, int b);

// Hairpin loop closed by pair (nuci, nucj); `size` = number of unpaired loop nucleotides.
// nuci1 = nt after i, nucj_1 = nt before j (interior mismatch). special_index >= 0 selects a
// tri/tetra/hexaloop bonus (see special_index()); -1 = generic.
int hairpin(int size, int nuci, int nuci1, int nucj_1, int nucj, int special_index);

// Two-loop (stack / bulge / internal) between outer pair (nuci,nucj) and inner pair (nucp,nucq).
// n1 = unpaired on the 5' side (between i and p), n2 = unpaired on the 3' side (between q and j).
int two_loop(int n1, int n2, int nuci, int nuci1, int nucj_1, int nucj,
             int nucp_1, int nucp, int nucq, int nucq1);

// Same as two_loop(), but with the outer/inner pair types precomputed by the caller (as
// pair_type(nuci,nucj) / pair_type(nucq,nucp)) instead of recomputed on every call. In
// two_loop_cost()'s nested enumeration, the outer type is invariant for an entire call and the
// inner type2 is invariant across the internal-loop's gap-option sub-loops, so callers that
// already hold these values can skip the redundant pair_type() work by calling this directly.
int two_loop_pt(int n1, int n2, int type, int nuci1, int nucj_1, int nucp_1, int type2, int nucq1);

// Multiloop terms (dangle model 0).
int multi_closing(int nuci, int nucj);   // closing pair (i,j) of a multiloop
int multi_branch(int nuci, int nucj);    // a stem (i,j) that is a branch of a multiloop
int multi_unpaired();                    // per unpaired base in a multiloop

// External loop terms (dangle model 0).
int external_paired(int nuci, int nucj); // a stem (i,j) in the external loop
int external_unpaired();                 // per unpaired base in the external loop

// Special hairpins: index of a tri(span 5)/tetra(span 6)/hexa(span 8) loop, else -1.
// `fullspan` is the closing pair + loop, read 5'->3' (e.g. "CUUCGG" for a tetraloop).
int special_index(const std::string& fullspan);

// Special hairpins of a given loop size (3=tri, 4=tetra, 6=hexa): list of
// (full span incl. closing pair, energy). Empty for other sizes. Used to match the
// ~22 known special loops against the design space without enumerating all loop sequences.
const std::vector<std::pair<std::string, int>>& special_hairpins(int size);

}  // namespace energy
}  // namespace ldclean

#endif  // LDCLEAN_ENERGY_H
