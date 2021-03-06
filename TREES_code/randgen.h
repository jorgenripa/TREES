//
// TREES -
// A TRait-based Eco-Evolutionary Simulation tool
// Copyright (C) 2017  Jörgen Ripa
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Contact: jorgen.ripa@biol.lu.se
//


#ifndef Species_rnd_h
#define Species_rnd_h

#include <vector>

#include "types.h"

double rand1(); // uniform [0,1)
double randn();
int poissrnd(double lambda);
int poissrndAppr(double lambda);
int binornd(int n, double p);
int binorndAppr(int n, double p);
int binoSmallP(int n, double p);
double doubleExp(double absmean);
seed_type randomize(); // also return the seed
void randseed(seed_type seed);
seed_type make_seed();
int weightedChoice( std::vector<double>& weights);
int weightedChoiceCumSum( std::vector<double>& cumWeights);

// class to generate a stream of random bits (true/false)
class bitGenerator {
    uint64_t bits;
    int bitsleft;
    void regenerate();

public:
    bitGenerator();
    inline bool nextBit(){
        if (bitsleft-- == 0) {
            regenerate();
        }
        bool reply = bits & 1;
        bits >>= 1;
        return reply;
    }
};

#endif
