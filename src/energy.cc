// energy.cc — see energy.h. Turner 2004 model, dangle model 0.
// Derived from the open-source ViennaRNA / LinearFold energy formulas; parameters in vienna/.
#include "energy.h"

#include <cmath>
#include <cstring>

// Turner 2004 parameters (open-source). energy_parameter.h #defines SPECIAL_HP; harmless here.
#ifdef SPECIAL_HP
#undef SPECIAL_HP
#endif
#include "vienna/energy_parameter.h"
#include "vienna/intl11.h"
#include "vienna/intl21.h"
#include "vienna/intl22.h"

namespace ldclean {
namespace energy {

namespace {
constexpr int MAXLOOP = 30;
// ViennaRNA pair types: 0:none 1:CG 2:GC 3:GU 4:UG 5:AU 6:UA.  Keeping one
// immutable table for both public queries avoids mutable initialization state and
// guarantees that pairable() and pair_type() cannot disagree.
constexpr unsigned char kPairType[5][5] = {
    {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 5},  // A-U
    {0, 0, 0, 1, 0},  // C-G
    {0, 0, 2, 0, 3},  // G-C / G-U
    {0, 6, 0, 4, 0},  // U-A / U-G
};

using HairpinList = std::vector<std::pair<std::string, int>>;

struct SpecialHairpinTables {
    HairpinList tri;
    HairpinList tet;
    HairpinList hex;
    HairpinList none;
};

inline int min2(int a, int b) { return a < b ? a : b; }
}  // namespace

void init() {}

bool pairable(int a, int b) {
    return (unsigned)a <= 4 && (unsigned)b <= 4 && kPairType[a][b] != 0;
}

// pair type (ViennaRNA): 0:none 1:CG 2:GC 3:GU 4:UG 5:AU 6:UA ; nucs A=1 C=2 G=3 U=4. Exposed
// (not just an internal helper) so hot callers can hoist it out of a loop and reuse the result
// across many two_loop_pt() calls instead of recomputing it every time (see energy.h).
int pair_type(int a, int b) {
    return ((unsigned)a <= 4 && (unsigned)b <= 4) ? kPairType[a][b] : 0;
}

int hairpin(int size, int nuci, int nuci1, int nucj_1, int nucj, int special_index) {
    int type = pair_type(nuci, nucj);
    int e = (size <= 30) ? hairpin37[size]
                         : (hairpin37[30] + (int)(lxc37 * std::log(size / 30.)));
    if (size < 3) return e;
    if (size == 4 && special_index > -1) return Tetraloop37[special_index];
    if (size == 6 && special_index > -1) return Hexaloop37[special_index];
    if (size == 3) {
        if (special_index > -1) return Triloop37[special_index];
        return e + (type > 2 ? TerminalAU37 : 0);
    }
    e += mismatchH37[type][nuci1][nucj_1];
    return e;
}

int two_loop_pt(int n1, int n2, int type, int nuci1, int nucj_1, int nucp_1, int type2, int nucq1) {
    int nl = n1 > n2 ? n1 : n2;
    int ns = n1 > n2 ? n2 : n1;

    if (nl == 0) return stack37[type][type2];  // stack

    if (ns == 0) {  // bulge
        int e = (nl <= MAXLOOP) ? bulge37[nl] : (bulge37[30] + (int)(lxc37 * std::log(nl / 30.)));
        if (nl == 1) {
            e += stack37[type][type2];
        } else {
            if (type > 2) e += TerminalAU37;
            if (type2 > 2) e += TerminalAU37;
        }
        return e;
    }

    if (ns == 1) {
        if (nl == 1) return int11_37[type][type2][nuci1][nucj_1];        // 1x1
        if (nl == 2) {                                                   // 2x1
            if (n1 == 1) return int21_37[type][type2][nuci1][nucq1][nucj_1];
            return int21_37[type2][type][nucq1][nuci1][nucp_1];
        }
        int e = (nl + 1 <= MAXLOOP) ? internal_loop37[nl + 1]
                                    : (internal_loop37[30] + (int)(lxc37 * std::log((nl + 1) / 30.)));
        e += min2(MAX_NINIO, (nl - ns) * ninio37);
        e += mismatch1nI37[type][nuci1][nucj_1] + mismatch1nI37[type2][nucq1][nucp_1];
        return e;
    }
    if (ns == 2 && nl == 2) return int22_37[type][type2][nuci1][nucp_1][nucq1][nucj_1];  // 2x2
    if (ns == 2 && nl == 3) {                                                            // 2x3
        int e = internal_loop37[5] + ninio37;
        e += mismatch23I37[type][nuci1][nucj_1] + mismatch23I37[type2][nucq1][nucp_1];
        return e;
    }
    // generic internal loop
    int u = nl + ns;
    int e = (u <= MAXLOOP) ? internal_loop37[u] : (internal_loop37[30] + (int)(lxc37 * std::log(u / 30.)));
    e += min2(MAX_NINIO, (nl - ns) * ninio37);
    e += mismatchI37[type][nuci1][nucj_1] + mismatchI37[type2][nucq1][nucp_1];
    return e;
}

int two_loop(int n1, int n2, int nuci, int nuci1, int nucj_1, int nucj,
             int nucp_1, int nucp, int nucq, int nucq1) {
    return two_loop_pt(n1, n2, pair_type(nuci, nucj), nuci1, nucj_1, nucp_1,
                        pair_type(nucq, nucp), nucq1);
}

// dangle model 0: no terminal-mismatch / dangle contributions in multi/external loops.
int multi_closing(int nuci, int nucj) {
    int type = pair_type(nucj, nuci);  // closing pair is reversed for the multiloop interior
    return (type > 2 ? TerminalAU37 : 0) + ML_intern37 + ML_closing37;
}
int multi_branch(int nuci, int nucj) {
    int type = pair_type(nuci, nucj);
    return (type > 2 ? TerminalAU37 : 0) + ML_intern37;
}
int multi_unpaired() { return ML_BASE37; }

int external_paired(int nuci, int nucj) {
    int type = pair_type(nuci, nucj);
    return (type > 2 ? TerminalAU37 : 0);
}
int external_unpaired() { return 0; }

const std::vector<std::pair<std::string, int>>& special_hairpins(int size) {
    // Since C++11, function-local static initialization is synchronized.  Publishing a
    // const aggregate after construction makes every subsequent concurrent lookup read-only.
    static const SpecialHairpinTables tables = [] {
        SpecialHairpinTables result;
        for (int i = 0; i < (int)(sizeof(Triloop37) / sizeof(int)); ++i)
            result.tri.emplace_back(std::string(Triloops + i * 6, 5),
                                    Triloop37[i]);  // 5-char span + space
        for (int i = 0; i < (int)(sizeof(Tetraloop37) / sizeof(int)); ++i)
            result.tet.emplace_back(std::string(Tetraloops + i * 7, 6),
                                    Tetraloop37[i]);  // 6-char span + space
        for (int i = 0; i < (int)(sizeof(Hexaloop37) / sizeof(int)); ++i)
            result.hex.emplace_back(std::string(Hexaloops + i * 9, 8),
                                    Hexaloop37[i]);  // 8-char span + space
        return result;
    }();
    if (size == 3) return tables.tri;
    if (size == 4) return tables.tet;
    if (size == 6) return tables.hex;
    return tables.none;
}

int special_index(const std::string& fullspan) {
    char* ts;
    if (fullspan.size() == 6) {  // tetraloop
        if ((ts = strstr(Tetraloops, fullspan.c_str()))) return (int)(ts - Tetraloops) / 7;
    } else if (fullspan.size() == 8) {  // hexaloop
        if ((ts = strstr(Hexaloops, fullspan.c_str()))) return (int)(ts - Hexaloops) / 9;
    } else if (fullspan.size() == 5) {  // triloop
        if ((ts = strstr(Triloops, fullspan.c_str()))) return (int)(ts - Triloops) / 6;
    }
    return -1;
}

}  // namespace energy
}  // namespace ldclean
