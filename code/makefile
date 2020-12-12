Compile=g++ -g -std=c++11
triangulation_MC: main.o triangulation.o metropolis.o measurement.o experiment.o PT.o
	${Compile} -o triangulation_MC main.o triangulation.o metropolis.o  PT.o measurement.o experiment.o
main.o: main.cpp triangulation.h
	${Compile} -o main.o -c main.cpp
PT.o: PT.cpp PT.h triangulation.h
	${Compile} -o PT.o -c PT.cpp
triangulation.o: triangulation.cpp triangulation.h
	${Compile} -o triangulation.o -c triangulation.cpp
metropolis.o: metropolis.cpp triangulation.h
	${Compile} -o metropolis.o -c metropolis.cpp
measurement.o: measurement.cpp triangulation.h
	${Compile} -o measurement.o -c measurement.cpp
experiment.o: experiment.cpp triangulation.h
	${Compile} -o experiment.o -c experiment.cpp
clean:
	rm *.o
