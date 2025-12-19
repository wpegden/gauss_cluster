CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pedantic

.PHONY: all clean

all: gauss_cluster

build: gauss_cluster

sa: gauss_cluster

# Build the simulator binary

gauss_cluster: src/main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f gauss_cluster
