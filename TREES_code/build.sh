#!/bin/sh

#  build.sh
#  TREES
#  Created by Jörgen Ripa

g++ -c -std=c++0x main.cpp
g++ -c -std=c++0x simulation.cpp
g++ -c -std=c++0x genetics.cpp
g++ -c -std=c++0x trait.cpp
g++ -c -std=c++0x space.cpp
g++ -c -std=c++0x fitness.cpp
g++ -c -std=c++0x mating.cpp
g++ -c -std=c++0x population.cpp
g++ -c -std=c++0x geneTracking.cpp
g++ -c -std=c++0x parameterFile.cpp
g++ -c -std=c++0x sample.cpp
g++ -c -std=c++0x simfile.cpp
g++ -c -std=c++0x fastmath.cpp
g++ -c -std=c++0x mytime.cpp
g++ -c -std=c++0x randgen.cpp
g++ -c -std=c++0x transform.cpp

g++ -o TREES main.o simulation.o genetics.o trait.o space.o fitness.o population.o mating.o geneTracking.o sample.o parameterFile.o simfile.o fastmath.o randgen.o mytime.o transform.o

rm *.o
