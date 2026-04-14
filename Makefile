CXX := clang++
CXXFLAGS := -std=c++20 -Wall -Wextra -Wpedantic -Werror -MMD -MP -Iinclude
OPTFLAGS := -O2 -g
PKG_CONFIG ?= pkg-config
CATCH2_CFLAGS := $(shell $(PKG_CONFIG) --cflags catch2)
CATCH2_LIBS := $(shell $(PKG_CONFIG) --libs catch2-with-main)

APP_SRCS := src/main.cpp src/instance.cpp src/preprocessing.cpp src/exact_solver.cpp src/rng.cpp src/rooted_tree.cpp src/decomposition/common.cpp src/decomposition/components.cpp src/decomposition/height_reduction.cpp src/decomposition/blocks.cpp src/decomposition/clusters.cpp src/decomposition/cells.cpp src/algorithms/one_point_five_approx.cpp
APP_OBJS := $(patsubst src/%.cpp,build/src/%.o,$(APP_SRCS))
APP_BIN := bin/tucvrp

TEST_SRCS := tests/instance_test.cpp tests/preprocessing_test.cpp tests/one_point_five_approx_test.cpp tests/rng_test.cpp tests/rooted_tree_test.cpp tests/decomposition_test.cpp
TEST_OBJS := $(patsubst tests/%.cpp,build/tests/%.o,$(TEST_SRCS)) $(patsubst src/%.cpp,build/src/%.o,$(filter-out src/main.cpp,$(APP_SRCS)))
TEST_BIN := bin/tests

APP_DEPS := $(patsubst %.o,%.d,$(APP_OBJS))
TEST_DEPS := $(patsubst %.o,%.d,$(TEST_OBJS))

.PHONY: all clean test run format

all: $(APP_BIN)

test: $(TEST_BIN)
	./$(TEST_BIN)

run: $(APP_BIN)
	./$(APP_BIN)

$(APP_BIN): $(APP_OBJS) | bin
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $^ -o $@

$(TEST_BIN): $(TEST_OBJS) | bin
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(CATCH2_CFLAGS) $^ $(CATCH2_LIBS) -o $@

build/src/algorithms/%.o: src/algorithms/%.cpp | build/src/algorithms
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) -c $< -o $@

build/src/decomposition/%.o: src/decomposition/%.cpp | build/src/decomposition
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) -c $< -o $@

build/src/%.o: src/%.cpp | build/src
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) -c $< -o $@

build/tests/%.o: tests/%.cpp | build/tests
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(CATCH2_CFLAGS) -c $< -o $@

build/src:
	mkdir -p $@

build/src/algorithms:
	mkdir -p $@

build/src/decomposition:
	mkdir -p $@

build/tests:
	mkdir -p $@

bin:
	mkdir -p $@

clean:
	rm -rf build bin

format:
	clang-format -i $(shell find include src tests -name '*.hpp' -o -name '*.cpp')

-include $(APP_DEPS) $(TEST_DEPS)
