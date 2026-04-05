# Galaxy N-body simulator (C++) - phase 1, single-threaded
# Requires C++11 or later
#
# Targets:
#   make / make all     -> galaxy_sim
#   make test           -> galaxy_tests + integration shell scripts
#   make check          -> alias for test (legacy)
#   make test_unit      -> doctest binary only

CXX = g++
# -MMD/-MP: track headers (e.g. config.hpp) so Config layout changes cannot leave stale .o files
# that disagree across TUs (undefined behavior: crashes in string ops during galaxy init).
CXXFLAGS = -std=c++11 -O2 -Wall -I. -Ivendor -MMD -MP

# Library objects (everything except main.cpp) — linked into galaxy_tests
LIB_SRCS = config.cpp compare_orchestration.cpp force_compare.cpp physics/registry.cpp physics/Newtonian/newtonian.cpp \
	physics/TPFCore/source_ansatz.cpp physics/TPFCore/field_evaluation.cpp \
	physics/TPFCore/derived_tpf_radial.cpp physics/TPFCore/provisional_readout.cpp \
	physics/TPFCore/tpf_core_package.cpp physics/TPFCore/v11_weak_field_correspondence.cpp \
	integrator.cpp galaxy_init.cpp init_conditions.cpp scenario_defaults.cpp resolved_scenario.cpp \
	git_provenance.cpp render_audit.cpp simulation.cpp output.cpp

LIB_OBJS = $(LIB_SRCS:.cpp=.o)

MAIN_SRC = main.cpp
MAIN_OBJ = main.o

SRCS = $(LIB_SRCS) $(MAIN_SRC)
OBJS = $(LIB_OBJS) $(MAIN_OBJ)

TARGET = galaxy_sim

# App-level + auto-discovered package unit/regression tests (integration stays shell-based).
# Note: 'tests/unit/*.cpp' does not match the pattern '*/tests/unit/*.cpp', so we combine wildcards + find.
APP_TEST_CPP := $(wildcard tests/unit/*.cpp tests/regression/*.cpp)
PKG_TEST_CPP := $(shell find physics -type f \( -path '*/tests/unit/*.cpp' -o -path '*/tests/regression/*.cpp' \) 2>/dev/null)
TEST_CPP := $(sort $(APP_TEST_CPP) $(PKG_TEST_CPP))
TEST_OBJS := $(TEST_CPP:.cpp=.o)

TEST_MAIN_SRC = tests/doctest_main.cpp
TEST_MAIN_OBJ = tests/doctest_main.o

TEST_TARGET = galaxy_tests

# Auto-generated dependency files (same stem as .o)
DEP_FILES := $(OBJS:.o=.d) $(TEST_OBJS:.o=.d) $(TEST_MAIN_OBJ:.o=.d)

.PHONY: all clean test check test_unit test_integration

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

$(TEST_TARGET): $(LIB_OBJS) $(TEST_OBJS) $(TEST_MAIN_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(LIB_OBJS) $(TEST_OBJS) $(TEST_MAIN_OBJ)

$(TEST_MAIN_OBJ): $(TEST_MAIN_SRC)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

test_unit: $(TARGET) $(TEST_TARGET)
	./$(TEST_TARGET)

test_integration: $(TARGET)
	./tests/run_integration.sh

# Full C++ suite: unit/regression binary + integration shell scripts (requires galaxy_sim).
test: $(TARGET) $(TEST_TARGET)
	./$(TEST_TARGET)
	./tests/run_integration.sh
check: test

clean:
	rm -f $(OBJS) $(TARGET) $(TEST_TARGET) $(TEST_OBJS) $(TEST_MAIN_OBJ) $(DEP_FILES)

# Must be last: .d files define rules like "config.o: ..." and must not become the default goal.
-include $(DEP_FILES)
