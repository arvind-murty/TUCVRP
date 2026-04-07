CXX := clang++
CXXFLAGS := -std=c++20 -Wall -Wextra -Wpedantic -Werror -MMD -MP -Iinclude
OPTFLAGS := -O2 -g
PKG_CONFIG ?= pkg-config
CATCH2_CFLAGS := $(shell $(PKG_CONFIG) --cflags catch2)
CATCH2_LIBS := $(shell $(PKG_CONFIG) --libs catch2-with-main)

APP_SRCS := src/main.cpp src/instance.cpp src/preprocessing.cpp src/solver.cpp
APP_OBJS := $(patsubst src/%.cpp,build/src/%.o,$(APP_SRCS))
APP_BIN := bin/tucvrp

TEST_SRCS := tests/instance_test.cpp tests/preprocessing_test.cpp
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

build/src/%.o: src/%.cpp | build/src
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) -c $< -o $@

build/tests/%.o: tests/%.cpp | build/tests
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(CATCH2_CFLAGS) -c $< -o $@

build/src:
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
