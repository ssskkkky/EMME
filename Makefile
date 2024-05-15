# Define the compiler
CXX = g++

# Define C++ compiler flags (feel free to customize)
CXXFLAGS = -Wall -std=c++17

# Define the main executable name
TARGET = emme

# Define all source files
SRCS = main.cpp Parameters.cpp solver.cpp

# Define all header files (usually only the main header)
HDRS = Parameters.h solver.h

# Build the executable
$(TARGET): $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS)

# Clean the project (removes the executable)
clean:
	rm -f $(TARGET)

# Default target to build the executable
.PHONY: all clean

all: $(TARGET)
