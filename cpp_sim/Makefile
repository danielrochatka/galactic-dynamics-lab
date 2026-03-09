# Galaxy N-body simulator (C++) - phase 1, single-threaded
# Requires C++11 or later

CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall -I.

SRCS = config.cpp physics/registry.cpp physics/Newtonian/newtonian.cpp physics/TPF/tpf_package.cpp integrator.cpp init_conditions.cpp simulation.cpp output.cpp main.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = galaxy_sim

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(TARGET)
