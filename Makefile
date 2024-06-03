# Define the compiler

CXX = g++

LAPACK_INCLUDE = $(shell pkg-config lapack --cflags)
BLASLAPCK_LIBS = $(shell pkg-config lapack --libs) $(shell pkg-config blas --libs)

# Define C++ compiler flags (feel free to customize)

DEBUG_FLAGS = -g -DEMME_DEBUG
OPT_FLAGS = -O3
CXXFLAGS = -Wall -std=c++20 -DEMME_EXPRESSION_TEMPLATE -DMULTI_THREAD -I./include

hash = $(subst ",\", $(shell ./build_info.sh hash))
date = $(subst ",\", $(shell ./build_info.sh time))

VPATH = src:include

ifdef DEBUG
CXXFLAGS +=$(DEBUG_FLAGS)
else
CXXFLAGS +=$(OPT_FLAGS)
endif

LD_FLAGS = $(BLASLAPCK_LIBS) -lgfortran

# Define the main executable name
TARGET = emme

# Define all source files
SRCS = $(notdir $(wildcard src/*.cpp))

OBJS = $(SRCS:.cpp=.o)

header_in_main = Grid.h JsonParser.h Matrix.h Parameters.h functions.h singularity_handler.h solver.h

all: $(TARGET)

Parameters.o: functions.h Timer.h
solver.o: Grid.h Matrix.h Parameters.h functions.h

# General Rules

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LD_FLAGS)

main.o: main.cpp $(header_in_main)
	$(CXX) $(CXXFLAGS) -c -o $@ $< -D$(hash) -D$(date)
$(filter-out main.o, $(OBJS)): %.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) $(LAPACK_INCLUDE) -c -o $@ $<

# Clean the project (removes the executable)
clean:
	rm -f $(TARGET) $(OBJS)

# Default target to build the executable
.PHONY: all clean remake test

test:
	make -C ./test/

remake:
	make clean;make -j
