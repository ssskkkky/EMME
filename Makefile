# Define the compiler
CXX = g++

# Define C++ compiler flags (feel free to customize)
CXXFLAGS = -Wall -std=c++17

# Define the main executable name
TARGET = emme

# Define all source files
SRCS = *.cpp

# Define all header files (usually only the main header)
HDRS = *.h

# Build the executable
$(TARGET): $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) -O2 -o $(TARGET) $(SRCS)

# Clean the project (removes the executable)
clean:
	rm -f $(TARGET)

# Default target to build the executable
.PHONY: all clean

all: $(TARGET)
