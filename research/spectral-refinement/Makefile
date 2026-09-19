CXX = g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic

build/diagmc: src/main.cpp src/diagmc.hpp src/conditional.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -static-libstdc++ -static-libgcc src/main.cpp -o $@

.PHONY: all test
all: build/diagmc

build/test_kernel: tests/kernel.cpp src/diagmc.hpp src/conditional.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -static-libstdc++ -static-libgcc tests/kernel.cpp -o $@

test: build/test_kernel
	./build/test_kernel
