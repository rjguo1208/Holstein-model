CXX = g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic

build/diagmc: src/main.cpp src/diagmc.hpp src/conditional.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -static-libstdc++ -static-libgcc src/main.cpp -o $@

.PHONY: all test
all: build/diagmc build/defect_diagmc build/defect_momentum_diagmc

build/defect_momentum_diagmc: src/defect_momentum_main.cpp src/defect_momentum.hpp src/defect.hpp src/conditional.hpp src/diagmc.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -static-libstdc++ -static-libgcc src/defect_momentum_main.cpp -o $@

build/test_defect_momentum: tests/defect_momentum_kernel.cpp src/defect_momentum.hpp src/defect.hpp src/conditional.hpp src/diagmc.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -static-libstdc++ -static-libgcc tests/defect_momentum_kernel.cpp -o $@

build/defect_diagmc: src/defect_main.cpp src/defect.hpp src/conditional.hpp src/diagmc.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -static-libstdc++ -static-libgcc src/defect_main.cpp -o $@

build/test_defect: tests/defect_kernel.cpp src/defect.hpp src/conditional.hpp src/diagmc.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -static-libstdc++ -static-libgcc tests/defect_kernel.cpp -o $@

build/test_kernel: tests/kernel.cpp src/diagmc.hpp src/conditional.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -static-libstdc++ -static-libgcc tests/kernel.cpp -o $@

test: build/test_kernel build/test_defect build/test_defect_momentum
	./build/test_kernel
	./build/test_defect
	./build/test_defect_momentum
