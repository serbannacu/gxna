CC = g++
CFLAGS = -g -std=c++17 -Wall -W -MMD -O3
LDFLAGS = -g

SOURCES = Args.cpp Distribution.cpp Statistics.cpp RandomPermutation.cpp GeneNetwork.cpp Experiment.cpp gxna.cpp
DEPENDS = $(addprefix bin/, $(SOURCES:.cpp=.d))

all: bin/gxna

clean:
	rm -rf bin
bin:
	mkdir -p bin

bin/gxna: $(addprefix bin/, $(SOURCES:.cpp=.o))

bin/%.o: %.cpp | bin
	$(CC) $(CFLAGS) -c $< -o $@

-include ${DEPENDS}
