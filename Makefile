CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++14 -Wall

Hydro-fluks: Hydro-fluks.cpp bj_flow.cpp bj_flow.h
	$(CXX) $(CXXFLAGS) -o $@ Hydro-fluks.cpp bj_flow.cpp

clean:
	rm -f Hydro-fluks Hydro-fluks-cn Hydro-fluks-rk4 tests/test_cn tests/test_rk4 compare/compare_solvers

.PHONY: clean

bj_flow_cn.o: bj_flow_cn.cpp bj_flow_cn.h bj_flow.h
	$(CXX) $(CXXFLAGS) -c bj_flow_cn.cpp -o $@

Hydro-fluks-cn: Hydro-fluks-cn.cpp bj_flow.cpp bj_flow_cn.cpp bj_flow_cn.h bj_flow.h
	$(CXX) $(CXXFLAGS) -o $@ Hydro-fluks-cn.cpp bj_flow.cpp bj_flow_cn.cpp

tests/test_cn: tests/test_cn.cpp bj_flow.cpp bj_flow_cn.cpp bj_flow_cn.h bj_flow.h
	$(CXX) $(CXXFLAGS) -o $@ tests/test_cn.cpp bj_flow.cpp bj_flow_cn.cpp

test: tests/test_cn
	cd tests && ./test_cn

bj_flow_rk4.o: bj_flow_rk4.cpp bj_flow_rk4.h bj_flow_cn.h bj_flow.h
	$(CXX) $(CXXFLAGS) -c bj_flow_rk4.cpp -o $@

Hydro-fluks-rk4: Hydro-fluks-cn.cpp bj_flow.cpp bj_flow_cn.cpp bj_flow_rk4.cpp bj_flow_rk4.h bj_flow_cn.h bj_flow.h
	$(CXX) $(CXXFLAGS) -DUSE_RK4 -o $@ Hydro-fluks-cn.cpp bj_flow.cpp bj_flow_cn.cpp bj_flow_rk4.cpp

tests/test_rk4: tests/test_rk4.cpp bj_flow.cpp bj_flow_cn.cpp bj_flow_rk4.cpp
	$(CXX) $(CXXFLAGS) -o $@ tests/test_rk4.cpp bj_flow.cpp bj_flow_cn.cpp bj_flow_rk4.cpp

test-rk4: tests/test_rk4
	cd tests && ./test_rk4

compare/compare_solvers: compare/compare_solvers.cpp bj_flow.cpp bj_flow_cn.cpp bj_flow_rk4.cpp
	$(CXX) $(CXXFLAGS) -o $@ compare/compare_solvers.cpp bj_flow.cpp bj_flow_cn.cpp bj_flow_rk4.cpp

# Green-function comparison of the three solvers: data in compare/results, figures in docs/figures.
# PREFIX selects the table set (default: the HRG tables built by tables/make_1n0_tables.py).
PREFIX ?= tables/1n0
compare: compare/compare_solvers
	mkdir -p compare/results docs/figures
	./compare/compare_solvers $(PREFIX) compare/results
	python3 compare/plot_green_comparison.py

# polynomial solver tables from the HRG trajectory tables of this folder
tables:
	cd tables && python3 make_1n0_tables.py

.PHONY: compare tables
