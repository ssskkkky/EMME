DEBUGFLAGS = -std=c++20 -ggdb -DEMME_DEBUG -I..
CXX = g++

# all tests
TESTS = test_json

all: $(TESTS)

# extra object dependences of each test
test_json: JsonParser.o

# extra header dependences of each .o file
test_json.o: ../JsonParser.h
JsonParser.o: ../JsonParser.h

# The rest should be seldom modified

$(TESTS): %: %.o
	$(CXX) -o $@ $+

%.o: %.cpp
	$(CXX) $(DEBUGFLAGS) -o $@ -c $<

%.cpp: ../%.cpp
	cp $< ./

.PHONY: all clean

clean:
	rm *.o test_json
