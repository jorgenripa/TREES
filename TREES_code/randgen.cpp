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


#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <chrono>
#include <cstdint>
#include <numeric>
using namespace std;

#include "randgen.h"

int myround(double x) { return (int)floor(x+0.5);}

unsigned randomize()
{
    // set seed as microseconds since Epoch.
    using namespace std::chrono;
    system_clock::time_point now = system_clock::now();
    
    unsigned micros = (unsigned)duration_cast<microseconds>(now.time_since_epoch()).count();

    srand(micros);
    std::cout << "Randomize seed : " << micros << '\n';
    return micros;
}

// unsigned long rand_count;

void randseed(unsigned seed) {
    srand(seed);
    std::cout << "Randseed set : " << seed << '\n';
//    rand_count = 0;
}

double rand1() // uniform [0,1)
{
    int r = rand();
    double C = int64_t(RAND_MAX)+1.0; // RAND_MAX is 2^31 in OS-X 10.9.5
//    ++rand_count;
    return r / C;//double(RAND_MAX+1);
}

/*unsigned long get_rand_count() {
    return rand_count;
}*/

double randn() {
    static bool iset = false;
    static double gset = 0;
    double fac, rsq, v1, v2;
    
    if (!iset) {
        do {
            v1 = 2.0*rand1() - 1.0;
            v2 = 2.0*rand1() - 1.0;
            rsq = v1*v1 + v2*v2;
        } while (rsq>=1.0 || rsq==0);
        fac = sqrt(-2.0*log(rsq)/rsq);
        gset = v1*fac;
        iset = true;
        return v2*fac;
    } else {
        iset = false;
        return gset;
    }
}

// A double exponential distributed random number,
// i.e. exponentially distributed in both + and - direction
// absmean is mean(abs(x))
double doubleExp(double absmean) {
    double r = rand1(); // r is recycled
    if (r<0.5) {
        // now r is uniform [0, 0.5)
        return -absmean*log(2*r);
    } else {
        // now r is uniform [0.5, 1)
        return absmean*log(2*r-1);
    }
}

int binoSmallP(int n, double p) {
    int r=0;
    double cumP=0;
    double lnPk = n*log(1-p); // log for precision, and avoiding zero Pk
    while( log(rand1()) > lnPk - log(1-cumP) ) { //% k is larger? "(1-rand)" instead of just "rand" to avoid error in some rand.m implementations
        r += 1;
        cumP += exp(lnPk);
        lnPk += log((double)((n-r+1))/r*p/(1-p));
    }
    return r;
}

int binornd(int n, double p)
{
    const double smallPlimit = 0.01;
    const double smallPlimitN = 100;
    int r;
    bool conjP;
    
    if (n==0 || p==0) {
        return 0;
    }
    conjP = p>0.5;
    if( conjP ) {
        p = 1-p;
    }
    if (p<=smallPlimit && n>smallPlimitN) {
        r = binoSmallP(n,p);
    } else {
        r = 0;
        for( int j=0; j<n; j++ ) {
            if( rand1() < p ) {
                r++;
            }
        }
    }
    if( conjP ) {
        r = n-r;
    }
    return r;
}

int poissrnd(double lambda)
{
    double emL, P, lnP;
    int r = 0;
    if (lambda==0) {
        return 0;
    }
    emL = exp(-lambda);
    P = rand1();
    if( emL>0 ){
        while( P>emL ) {
            P = P*rand1();
            r++;
        }
    } else {
        lnP = -log(P);
        while( lnP<lambda ) {
            lnP = lnP-log(rand1());
            r++;
        }
    }
    return r;
}


int binorndAppr(int n, double p)
{
    int nLimit = 10; // below this limit, don't bother approximations, use exact method
    double NormalLimit = 150;
    double PoissonLimit = 0.035; // These limits assure the total misplaced prob mass is below 1%
    
    double npq = n*p*(1-p);
    if( npq > NormalLimit ){
        return max(0, myround(n*p + sqrt(npq)*randn()));
    } else if( n>nLimit && p<PoissonLimit ) {
        return min(n, poissrnd(n*p));
    } else {
        return binornd(n,p);
    }
}

int poissrndAppr(double lambda) {
    double large_limit = 160; // A limit of 160 means less than 1% of the probability mass will be misplaced
    int x;
    if( lambda>large_limit ){
        x = std::max(0, myround(lambda + sqrt(lambda)*randn()));
    }
    else {
        x = poissrnd(lambda);
    }
    return x;
}

vector<double> cumWeights;

int weightedChoice( vector<double>& weights) {
    cumWeights.resize(weights.size());
    // calculate cumulative sum:
    std::partial_sum(weights.begin(), weights.end(), cumWeights.begin());
    return weightedChoiceCumSum( cumWeights );
}

int weightedChoiceCumSum( vector<double>& cumW ) {
    // choose random item (weighted by weigths)
    double z = rand1()*cumW.back();
    vector<double>::iterator cwp = cumW.begin();
    while (z>*cwp) {
        ++cwp;
    }
    return (int) (cwp - cumW.begin());
}

// This is inline for speed:
/* void bitGenerator::regenerate() {
    bits = (unsigned)rand1()*65536;
    bitsleft = 16;
}*/

bitGenerator::bitGenerator(){
    regenerate();
}

void bitGenerator::regenerate() {
    bits = (unsigned)(rand1()*65536);
    bitsleft = 16;
}
bool bitGenerator::nextBit() {
    //return rand1()<0.5;
    if (bitsleft==0) {
        regenerate();
    }
    bool reply = bits & 1;
    bits >>= 1;
    --bitsleft;
    return reply;
}


// This is inline for speed:
/* bool bitGenerator::nextBit() {
    if (bitsleft==0) {
        regenerate();
    }
    bool reply = bits % 2;
    bits >>= 1;
    --bitsleft;
    return reply;
} */

