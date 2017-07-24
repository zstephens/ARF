CXX=g++
CXXARGS=-Wall -O3 -g
all: 01_jf2seeds 05_extendExactMatches

01_jf2seeds: 01_jf2seeds.o
	$(CXX) $(CXXARGS) -o 01_jf2seeds 01_jf2seeds.o

01_jf2seeds.o: 01_jf2seeds.cpp
	$(CXX) $(CXXARGS) -c -o 01_jf2seeds.o 01_jf2seeds.cpp

05_extendExactMatches: 05_extendExactMatches.o
	$(CXX) $(CXXARGS) -o 05_extendExactMatches 05_extendExactMatches.o

05_extendExactMatches.o: 05_extendExactMatches.cpp
	$(CXX) $(CXXARGS) -c -o 05_extendExactMatches.o 05_extendExactMatches.cpp

clean:
	rm 01_jf2seeds
	rm 01_jf2seeds.o
	rm 05_extendExactMatches
	rm 05_extendExactMatches.o