# Define the compiler

# DEBUGFLAGS = -g -DEMME_DEBUG 

CXX = g++

LAPACK_INCLUDE = $(shell pkg-config lapack --cflags)
BLASLAPCK_LIBS = $(shell pkg-config lapack --libs) $(shell pkg-config blas --libs)

# Define C++ compiler flags (feel free to customize)

ifdef DEBUG
CXXFLAGS = -Wall -std=c++20 -DEMME_EXPRESSION_TEMPLATE -DMULTI_THREAD -g -DEMME_DEBUG
else
CXXFLAGS = -Wall -std=c++20 -DEMME_EXPRESSION_TEMPLATE -DMULTI_THREAD -O3
endif


LD_FLAGS = $(BLASLAPCK_LIBS) -lgfortran

# Define the main executable name
TARGET = emme

# Define all source files
SRCS = $(wildcard *.cpp)

OBJS = $(SRCS:.cpp=.o)

# Define all header files (usually only the main header)
HDRS = *.h


# Build the executable

$(TARGET): $(OBJS) $(HDRS)
	$(CXX) -o $@ $(OBJS) $(STATIC_LIBS) $(LD_FLAGS) $(DEBUGFLAGS) 


%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS)  $(LAPACK_INCLUDE) -c $< -o $@

# Clean the project (removes the executable)
clean:
	rm -f $(TARGET) $(OBJS)

# Default target to build the executable
.PHONY: all clean remake test

all: $(TARGET)

test:
	make -C ./test/

remake:
	make clean;make -j
