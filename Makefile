CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++11 -Wall

Hydro-fluks: Hydro-fluks.cpp bj_flow.cpp bj_flow.h
	$(CXX) $(CXXFLAGS) -o $@ Hydro-fluks.cpp bj_flow.cpp

clean:
	rm -f Hydro-fluks

.PHONY: clean

bj_flow_cn.o: bj_flow_cn.cpp bj_flow_cn.h bj_flow.h
	$(CXX) $(CXXFLAGS) -c bj_flow_cn.cpp -o $@

Hydro-fluks-cn: Hydro-fluks-cn.cpp bj_flow.cpp bj_flow_cn.cpp bj_flow_cn.h bj_flow.h
	$(CXX) $(CXXFLAGS) -o $@ Hydro-fluks-cn.cpp bj_flow.cpp bj_flow_cn.cpp

tests/test_cn: tests/test_cn.cpp bj_flow.cpp bj_flow_cn.cpp bj_flow_cn.h bj_flow.h
	$(CXX) $(CXXFLAGS) -o $@ tests/test_cn.cpp bj_flow.cpp bj_flow_cn.cpp

test: tests/test_cn
	cd tests && ./test_cn
