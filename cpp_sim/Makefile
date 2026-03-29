# Galaxy N-body simulator (C++) - phase 1, single-threaded
# Requires C++11 or later

CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall -I.

SRCS = config.cpp force_compare.cpp physics/registry.cpp physics/Newtonian/newtonian.cpp physics/TPFCore/source_ansatz.cpp physics/TPFCore/field_evaluation.cpp physics/TPFCore/derived_tpf_radial.cpp physics/TPFCore/provisional_readout.cpp physics/TPFCore/tpf_core_package.cpp integrator.cpp galaxy_init.cpp init_conditions.cpp render_audit.cpp simulation.cpp output.cpp main.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = galaxy_sim

TEST_TPF_OBJS = test_tpf_regression.o physics/TPFCore/source_ansatz.o physics/TPFCore/derived_tpf_radial.o

.PHONY: all clean test_tpf_regression check

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

test_tpf_regression: $(TEST_TPF_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(TEST_TPF_OBJS)

check: test_tpf_regression
	./test_tpf_regression

clean:
	rm -f $(OBJS) $(TARGET) test_tpf_regression $(TEST_TPF_OBJS)
