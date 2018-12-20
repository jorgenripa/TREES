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


#include <math.h>
#include <iostream>
#include "fastmath.h"

const int Fastexp::MAXX = 400;
double table[Fastexp::MAXX*10+1]; // double is actually faster than float on modern platforms
double difftable[Fastexp::MAXX*10];
// Fill a table of values from exp(-400) to exp(0)
void fillExpTable() {
    for (int i=0; i<=Fastexp::MAXX*10; ++i) {
        table[i] = exp(-i/10.0);
        if (i>0) {
            difftable[i-1]=table[i]-table[i-1];
        }
    }
}

Fastexp::Fastexp() {
    etable = table;
    edifftable = difftable;
}

