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


#ifndef __Species__simfile__
#define __Species__simfile__

#include <stdio.h>
#include <fstream>
#include <cstdint>

#include "geneTracking.h"

// class for saving data in binary format
// protected inheritance to prevent other output than
// of supported types
class Simfile : protected std::ofstream {
public:
    Simfile(std::string filename);
    ~Simfile();
    void write(float x);
    void write(int i);
    void write(int64_t i);
    void write(idType id);
    void writeArray( float* A, int size);
    void writeArray( double* A, int size);
    void writeArray( int* A, int size);
    void writeArray( idType* A, int size);
    void writeString(std::string s);
};

#endif /* defined(__Species__simfile__) */
