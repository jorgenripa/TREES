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

#include <cassert>

#include "bit_vector.h"

bit_vector::bit_vector() { // default constructor
    current_int = 0;
    next_position = 0;
}

void bit_vector::reserve(size_type bit_count) {
    size_type int_count = bit_count/(sizeof(store_type)*8) + 1;
    std::vector<store_type>:: reserve(int_count);
}

void bit_vector::push_back(bool bit) {
    if (bit) {
        store_type mask = store_type(1)<<next_position;
        current_int |= mask;
    }
    next_position++;
    if (next_position==sizeof(store_type)*8) {
        std::vector<store_type>::push_back(current_int); // store in base class
        current_int = 0;
        next_position = 0;
    }
}

bit_vector::size_type bit_vector::get_total_bit_count() {
    return size()*sizeof(store_type)*8 + next_position;
}

bit_vector::size_type bit_vector::get_total_byte_count() {
    return size()*sizeof(store_type);
}

void bit_vector::write_to_file(oSimfile& osf) {
    // Write size of bitstream, followed by the bits themselves
    osf.write<uint32_t>((uint32_t)get_total_bit_count());
    osf.write<uint32_t>((uint32_t)sizeof(store_type));
    if(size()>0) {
        osf.writeArray<store_type>(&at(0), (int)size());
    }
    if (next_position>0) {
        osf.write<store_type>(current_int);
    }
}

void bit_vector::read_from_file(iSimfile &isf) {
    clear();
    uint32_t total_bits = isf.read<uint32_t>();
    uint32_t chunk_size = isf.read<uint32_t>(); // size of sore_type (in bytes)
    assert(chunk_size==sizeof(store_type));
    uint32_t total_bytes = total_bits / (sizeof(store_type)*8);
    next_position = total_bits % (sizeof(store_type)*8); // rest bits
    if(total_bytes>0) {
        isf.readArray<store_type>(*this, total_bytes);
    }
    if (next_position>0) {
        current_int = isf.read<store_type>();
    } else {
        current_int = 0;
    }
}

bool bit_vector::get_bit(size_type bit_pos) {
    size_type int_pos = bit_pos/(sizeof(store_type)*8);
    store_type the_int;
    if (int_pos<size()) {
        the_int = at(int_pos);
    } else {
        the_int = current_int;
    }
    int bit = bit_pos % (sizeof(store_type)*8);
    return (the_int & store_type(1)<<bit)>0;
}
