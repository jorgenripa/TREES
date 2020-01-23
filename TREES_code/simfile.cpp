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
#include <sstream>
#include <assert.h>
#include "simfile.h"

void endian_test() {
    // Test for Little Endian:
    unsigned char SwapTest[2] = { 1, 0 };
    short w = *(short *) SwapTest;
    assert( w == 1 );
}

oSimfile::oSimfile(std::string filename) {
    endian_test();
    open(filename.c_str(), std::fstream::binary | std::fstream::out);
    imbue(std::locale::classic());
    write<uint32_t>(current_file_version);
}

void oSimfile::close() {
    std::ofstream::close();
}

oSimfile::~oSimfile() {
    close();
}

void oSimfile::writeString(const std::string& s) {
    std::ofstream::write(s.c_str(), s.length());
    put(0);
}

iSimfile::iSimfile(std::string filename) {
    endian_test();
    open(filename.c_str(), std::fstream::binary | std::fstream::in);
    if (!is_open()) {
        std::cout << "Failed to open file " << filename << std::endl;
        exit(1);
    }
    imbue(std::locale::classic());
    int file_version = read<uint32_t>();
    if (file_version!=current_file_version ) {
        std::cout << "File error. File " << filename << " is version " << file_version << " instead of " << current_file_version << std::endl;
        exit(1);
    }
}

iSimfile::~iSimfile() {
    close();
}

std::string iSimfile::readString() {
    std::stringstream buffer;
    std::ifstream::get(*buffer.rdbuf(), char(0));
    // get the null character:
    std::ifstream::get();
    return buffer.str();
}



/*
void Simfile::write(float x) {
    std::ofstream::write((char*)&x, sizeof(float));
}

void Simfile::write(int i) {
    std::ofstream::write((char*)&i, sizeof(int));
}

void Simfile::write(int64_t i) {
    std::ofstream::write((char*)&i, sizeof(int64_t));
}

void Simfile::write(id_type id) {
    std::ofstream::write((char*)&id, sizeof(id_type));
}
*/


/*
void Simfile::writeArray( signed char* A, int size) {
    std::ofstream::write((char*)A, sizeof(signed char)*size);
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

void Simfile::writeArray( id_type* A, int size) {
    std::ofstream::write((char*)A, sizeof(id_type)*size);
}
*/
