CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++11 -Wall

Hydro-fluks: Hydro-fluks.cpp bj_flow.cpp bj_flow.h
	$(CXX) $(CXXFLAGS) -o $@ Hydro-fluks.cpp bj_flow.cpp

clean:
	rm -f Hydro-fluks

.PHONY: clean
