#!/bin/sh

#  build.sh
#  TREES
#  Created by Jörgen Ripa on 26/6 -14.

g++ -c -std=c++0x main.cpp
g++ -c -std=c++0x simulation.cpp
g++ -c -std=c++0x population.cpp
g++ -c -std=c++0x trait.cpp
g++ -c -std=c++0x transform.cpp
g++ -c -std=c++0x fitness.cpp
g++ -c -std=c++0x mating.cpp
g++ -c -std=c++0x space.cpp
g++ -c -std=c++0x genetics.cpp
g++ -c -std=c++0x geneTracking.cpp
g++ -c -std=c++0x genotype_phenotype_map.cpp
g++ -c -std=c++0x parameterFile.cpp
g++ -c -std=c++0x simfile.cpp
g++ -c -std=c++0x bit_vector.cpp
g++ -c -std=c++0x fastmath.cpp
g++ -c -std=c++0x randgen.cpp
g++ -c -std=c++0x mytime.cpp

g++ -o TREES main.o simulation.o population.o trait.o transform.o fitness.o mating.o space.o genetics.o geneTracking.o genotype_phenotype_map.o parameterFile.o simfile.o bit_vector.o fastmath.o randgen.o mytime.o

rm *.o
