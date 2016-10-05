CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11 -Ofast -fopenmp -msse4.2
LIBS =

OBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))

main.out: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o main.out $(LIBS)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LIBS)

profile: $(OBJS)
	$(CXX) $(CXXFLAGS) -pg $(OBJS) -o main.out $(LIBS)

.PHONY: clean
clean:
	\rm *.o *.out profile
