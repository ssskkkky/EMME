DEBUGFLAGS = -Wall -std=c++20 -ggdb -DEMME_DEBUG -I../include
CXX = /lustre/opt/latest-compiler/bin/g++
# CXX = g++

# all tests
TESTS = test_json test_integrator

lib_include_path = $(shell realpath .)/../include
lib_source_path = $(shell realpath .)/../src

all: $(TESTS)

# extra object dependences of each test
test_json: JsonParser.o

# extra header dependences of each .o file
test_json.o: $(lib_include_path)/JsonParser.h
JsonParser.o: $(lib_include_path)/JsonParser.h
test_integrator.o: $(lib_include_path)/solver_pic.h

# The rest should be seldom modified

$(TESTS): %: %.o
	$(CXX) -o $@ $+

%.o: %.cpp
	$(CXX) $(DEBUGFLAGS) -o $@ -c $<

%.cpp: $(lib_source_path)/%.cpp
	cp $< ./

.PHONY: all clean

clean:
	rm -f *.o $(TESTS)
