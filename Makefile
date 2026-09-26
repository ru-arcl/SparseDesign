# SparseDesign — standalone source build.
CXX ?= g++
CXXFLAGS ?= -std=c++11 -O3 -flto -Wall -fopenmp
# The square arena is the speed-oriented production default.  Omitting this flag selects the
# forward-only packed arena, exposed below as an explicitly named low-memory build.
DEFAULT_LAYOUT_FLAGS ?= -DLDCLEAN_SQUARE_MEMOS
SRC = src/main.cc src/fold_turner.cc src/energy.cc src/dfa.cc src/codon_table.cc
HDR = src/codon_table.h src/dfa.h src/energy.h src/fold_turner.h $(wildcard src/vienna/*.h)
BIN = sparsedesign
LOWMEM_BIN = sparsedesign-lowmem-packed
LEGACY_BIN = lineardesign-clean
LEGACY_LOWMEM_BIN = lineardesign-clean-lowmem-packed

all: $(BIN) $(LEGACY_BIN)

$(BIN): $(SRC) $(HDR)
	$(CXX) $(CXXFLAGS) $(DEFAULT_LAYOUT_FLAGS) $(SRC) -o $(BIN)

$(LOWMEM_BIN): $(SRC) $(HDR)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(LOWMEM_BIN)

# Preserve executable names used by the archived research scripts.
$(LEGACY_BIN): $(BIN)
	cp $(BIN) $(LEGACY_BIN)

$(LEGACY_LOWMEM_BIN): $(LOWMEM_BIN)
	cp $(LOWMEM_BIN) $(LEGACY_LOWMEM_BIN)

.PHONY: all lowmem-packed test test-layouts test-asan test-ubsan test-sanitize \
	test-lowmem-asan test-lowmem-ubsan clean

lowmem-packed: $(LOWMEM_BIN) $(LEGACY_LOWMEM_BIN)

test: $(BIN)
	@$(CXX) $(CXXFLAGS) $(DEFAULT_LAYOUT_FLAGS) test/test_codon.cc src/codon_table.cc -o /tmp/ldc_tc && (cd . && /tmp/ldc_tc)
	@$(CXX) $(CXXFLAGS) $(DEFAULT_LAYOUT_FLAGS) test/test_dfa.cc src/dfa.cc src/codon_table.cc -o /tmp/ldc_td && (cd . && /tmp/ldc_td)
	@$(CXX) $(CXXFLAGS) $(DEFAULT_LAYOUT_FLAGS) test/test_fold_simple.cc src/fold_simple.cc src/dfa.cc src/codon_table.cc -o /tmp/ldc_tf && (cd . && /tmp/ldc_tf)
	@$(CXX) $(CXXFLAGS) $(DEFAULT_LAYOUT_FLAGS) -pthread test/test_energy_threadsafe.cc src/energy.cc -o /tmp/ldc_tet && (cd . && /tmp/ldc_tet)
	@$(CXX) $(CXXFLAGS) $(DEFAULT_LAYOUT_FLAGS) test/test_sparse_exact.cc src/fold_turner.cc src/energy.cc src/dfa.cc src/codon_table.cc -o /tmp/ldc_tse && (cd . && /tmp/ldc_tse)
	@$(CXX) $(CXXFLAGS) $(DEFAULT_LAYOUT_FLAGS) test/test_motif_dfa.cc src/fold_turner.cc src/energy.cc src/dfa.cc src/codon_table.cc -o /tmp/ldc_tmd && (cd . && /tmp/ldc_tmd)
	@$(CXX) $(CXXFLAGS) $(DEFAULT_LAYOUT_FLAGS) test/test_turner_multiloop.cc src/fold_turner.cc src/energy.cc src/dfa.cc src/codon_table.cc -o /tmp/ldc_tml && (cd . && /tmp/ldc_tml)
	@$(CXX) $(CXXFLAGS) $(DEFAULT_LAYOUT_FLAGS) test/test_dense_wavefront.cc src/fold_turner.cc src/energy.cc src/dfa.cc src/codon_table.cc -o /tmp/ldc_tdw && (cd . && /tmp/ldc_tdw)
	@$(CXX) $(CXXFLAGS) $(DEFAULT_LAYOUT_FLAGS) test/test_sparse_kernel.cc -o /tmp/ldc_tsk && /tmp/ldc_tsk
	@$(CXX) $(CXXFLAGS) $(DEFAULT_LAYOUT_FLAGS) test/test_right_normal_cells.cc src/energy.cc src/dfa.cc src/codon_table.cc -o /tmp/ldc_trnc && /tmp/ldc_trnc
	@bash test/test_cli.sh ./$(BIN)
	@python3 test/test_infer_lambda.py ./$(BIN)

test-layouts:
	@bash test/test_layouts.sh

test-asan:
	@bash test/test_sanitizers.sh asan $(DEFAULT_LAYOUT_FLAGS)

test-ubsan:
	@bash test/test_sanitizers.sh ubsan $(DEFAULT_LAYOUT_FLAGS)

test-sanitize: test-asan test-ubsan

test-lowmem-asan:
	@bash test/test_sanitizers.sh asan

test-lowmem-ubsan:
	@bash test/test_sanitizers.sh ubsan

clean:
	rm -f $(BIN) $(LOWMEM_BIN) $(LEGACY_BIN) $(LEGACY_LOWMEM_BIN)
