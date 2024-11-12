# The recommended build uses cmake.
# This file is provided as a workaround in case cmake is not available.
# It should work on MacOS and Linux.
# It does not build the unit test binaries.

CC = g++
CFLAGS = -std=c++17 -Wall -W -MMD -O3
LDFLAGS =

SOURCES = Args.cpp Distribution.cpp Experiment.cpp GeneNetwork.cpp RandomPermutation.cpp Statistics.cpp gxna.cpp
DEPENDS = $(addprefix build/, $(SOURCES:.cpp=.d))
OBJS = $(addprefix build/, $(SOURCES:.cpp=.o))

all: build/gxna

clean:
	rm -rf build
build:
	mkdir -p build

build/gxna: $(OBJS)

build/%.o: */%.cpp | build
	$(CC) $(CFLAGS) -Iinclude -c $< -o $@

-include $(DEPENDS)
