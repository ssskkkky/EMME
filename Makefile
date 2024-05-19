# Define the compiler
CXX = g++

# Define C++ compiler flags (feel free to customize)
CXXFLAGS = -Wall -std=c++20 -O2

# LD_FLAGS = -L./deps/lapack/lib -llapack -lblas -Wl,-v
LD_FLAGS =  $(shell pkg-config lapack --libs) $(shell pkg-config blas --libs) -lgfortran
# STATIC_LIBS = ./deps/lapack/lib/liblapack.a ./deps/lapack/lib/libblas.a 

# Define the main executable name
TARGET = emme

# Define all source files
SRCS = $(wildcard *.cpp)

OBJS = $(SRCS:.cpp=.o)

# Define all header files (usually only the main header)
HDRS = *.h

# Build the executable
$(TARGET): $(OBJS) $(HDRS)
	$(CXX) -o $@ $(OBJS) $(STATIC_LIBS) $(LD_FLAGS) 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(shell pkg-config lapack --cflags) -c $< -o $@

# Clean the project (removes the executable)
clean:
	rm -f $(TARGET) $(OBJS)

# Default target to build the executable
.PHONY: all clean remake

all: $(TARGET)

remake:
	make clean;make -j
