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


#include <iostream>
#include <assert.h>
#include "simfile.h"

Simfile::Simfile(std::string filename) {
    // Test for Little Endian:
    unsigned char SwapTest[2] = { 1, 0 };
    short w = *(short *) SwapTest;
    assert( w == 1 );
        
    open(filename.c_str(), std::fstream::binary | std::fstream::out);
    imbue(std::locale::classic());
}

Simfile::~Simfile() {
    close();
}

void Simfile::write(float x) {
    std::ofstream::write((char*)&x, sizeof(float));
}

void Simfile::write(int i) {
    std::ofstream::write((char*)&i, sizeof(int));
}

void Simfile::write(int64_t i) {
    std::ofstream::write((char*)&i, sizeof(int64_t));
}

void Simfile::write(idType id) {
    std::ofstream::write((char*)&id, sizeof(idType));
}

void Simfile::writeArray( float* A, int size) {
    std::ofstream::write((char*)A, sizeof(float)*size);
}

void Simfile::writeArray( double* A, int size) {
    // Convert all doubles to float before writing:
    float* x = new float[size];
    for (int i=0; i<size; ++i) {
        x[i] = (float)A[i];
    }
    writeArray(x, size);
    delete [] x;
}

void Simfile::writeArray( int* A, int size) {
    std::ofstream::write((char*)A, sizeof(int)*size);
}

void Simfile::writeArray( idType* A, int size) {
    std::ofstream::write((char*)A, sizeof(idType)*size);
}

void Simfile::writeString(std::string s) {
    std::ofstream::write(s.c_str(), s.length());
    put('\n');
    put(0);
}
