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

#ifndef bit_vector_hpp
#define bit_vector_hpp

#include <stdio.h>
#include <vector>
#include <cstdint>
#include "simfile.h"

typedef uint64_t store_type;

class bit_vector : protected std::vector<store_type> {
    protected:
    store_type current_int;
    int next_position;
    
    public:
    bit_vector(); // default constructor
    void reserve(size_type bit_count);
    void push_back(bool bit);
    size_type get_total_bit_count();
    size_type get_total_byte_count();
    void write_to_file(oSimfile& osf);
    void read_from_file(iSimfile& isf);
    bool get_bit(size_type i);
};
#endif /* bit_vector_hpp */
