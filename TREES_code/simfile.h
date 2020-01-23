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
#include <iostream>
#include <fstream>
#include <cstdint>

#include "geneTracking.h"

static const uint32_t current_file_version = 4; // as of TREES 1.3

// class for saving data in binary format
// protected inheritance to prevent non-supported output
class oSimfile : protected std::ofstream {
public:
    oSimfile(std::string filename);
    ~oSimfile();
    template<typename WriteType> void write(const WriteType data) {
        std::ofstream::write((char*)&data, sizeof(WriteType));

    }
    template<typename WriteType>void writeArray( const WriteType * A, int size) {
        std::ofstream::write((char*)A, sizeof(WriteType)*size);
    }
    template<typename WriteType>void writeCArray( const WriteType * A, int size) {
        std::ofstream::write((char*)A, sizeof(WriteType)*size);
    }
    // This also stores vector length:
    template<typename WriteType>void writeVector( const std::vector<WriteType>& v) {
        write<size_type>((size_type)v.size());
        std::ofstream::write((char*)&v[0], sizeof(WriteType)*v.size());
    }
    
    void writeString(const std::string& s);
    
    void close();
    
/*    void write(float x);
    void write(int i);
    void write(int64_t i);
    void write(id_type id);
 */
/*
    void writeArray( signed char* A, int size);
    void writeArray( int* A, int size);
    void writeArray( float* A, int size);
    void writeArray( double* A, int size);
    void writeArray( id_type* A, int size);
 */
};

/*
template<typename T> oSimfile& operator <<(oSimfile& osf, std::vector<T>& v) {
    osf.write<int>((int)v.size());
    for (int i=0; i<v.size(); ++i) {
        osf << v[i]; // only works for some supported types
    }
    return osf;
}
*/

class iSimfile : protected std::ifstream {
public:
    iSimfile(std::string filename);
    ~iSimfile();
//    template<typename ReadType> void read(ReadType& data) {
//        std::ifstream::read((char*)&data, sizeof(ReadType));
//
//    }
    template<typename ReadType> ReadType read() {
        ReadType data;
        std::ifstream::read((char*)&data, sizeof(ReadType));
        return data;
    }
    template<typename ReadType>void readArray( std::vector<ReadType>& A, int size) {
        A.assign(size,(ReadType)0);
        std::ifstream::read((char*)&A[0], sizeof(ReadType)*size);
    }
    template<typename ReadType>void readCArray( ReadType* Ap, int count) {
        std::ifstream::read((char*)Ap, sizeof(ReadType)*count);
    }
    template<typename ReadType>void readVector(std::vector<ReadType>& v) {
        size_type vlength = read<size_type>();
        readArray<ReadType>(v, vlength);
    }
    std::string readString();
};

template<typename T> iSimfile& operator >>(iSimfile& isf, std::vector<T>& v) {
    int count = isf.read<int>();
    v.assign(count,T());
    for (int i=0; i<count; ++i) {
        isf >> v[i]; // only works for some supported types
    }
}

#endif /* defined(__Species__simfile__) */
