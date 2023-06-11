CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11

main.exe: main.o random.o
	$(CXX) random.o main.o -o main.exe

main.o: main.cpp random.h TSP.h
	$(CXX) -c main.cpp -o main.o $(CXXFLAGS)

random.o: random.cpp random.h
	$(CXX) -c random.cpp -o random.o $(CXXFLAGS)

clean:
	rm *.o main.exe seed.out

