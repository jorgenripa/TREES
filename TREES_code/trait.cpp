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

#include <cstdint>
#include <cassert>

#include "trait.h"
#include "population.h"
#include "genotype_phenotype_map.h"
#include "randgen.h"

// protected constructor, used by subclasses:
Trait::Trait(std::string& tname, Population& p) :
Matrix<traitType>(0,0), name(tname), pop(p), genes(p.getGenetics()) {
    GP_map = NULL;
    Xinit = 0;
}

Trait::Trait(std::string& tname, Population& p, ParameterFile& pf) :
Matrix<traitType>(0,0), name(tname), pop(p), genes(p.getGenetics()) {
    
    set_dims(pf.getPositiveInt("dimensions"));
    int lociPerDim = pf.getPositiveInt("loci_per_dim");
    Xinit =pf.getDouble("initial_value");
    GP_map = genes.create_GP_map(*this,lociPerDim, pf);
    
    // Reserve memory space for phenotypes (this will be expanded when necessary):
    reserve(get_dims(),10000);
}

bool Trait::is_constant() { return false; }

void Trait::initialize() {
    assign(get_dims(), pop.size(), 0);
    GP_map->initialize(Xinit); // this initializes the genes, not trait values
}

Trait::~Trait() {
    if (GP_map) {
        delete GP_map;
    }
}

void Trait::generatePhenotypes() {
    GP_map->generate_phenotypes();
}

traitType& Trait::traitValue(int individual, int dim) {
    assert(dim<get_dims() && individual<get_N());
    return (*this)(dim,individual);
}

/*
void Trait::add_to_checkpoint(Checkpoint &cp) {
    // Checkpoints don't store trait values.
    // But, they may contain GP_map parameters:
    GP_map->add_to_checkpoint(cp);
}

void Trait::read_to_checkpoint(Checkpoint &cp, iSimfile &isf) {
    GP_map->read_to_checkpoint(cp,isf);
}

int Trait::resume_from_checkpoint(Checkpoint &cp, int dataIndex) {
    dataIndex = GP_map->resume_from_checkpoint(cp,dataIndex);
    return dataIndex;
}
*/

//////////////////////////////////////////////////////////
// TraitConstant
///////////////////////////////////////////////////

TraitConstant::TraitConstant(std::string& name, Population& p, ParameterFile& pf):
Trait(name,p) {
    
    // A constant trait is implemented as a trait with zero loci:
    set_dims(pf.getPositiveInt("dimensions"));
    Xinit =pf.getDouble("initial_value");
    
    // Reserve memory space for phenotypes (this will be expanded when necessary):
    reserve(get_dims(),10000);
}

bool TraitConstant::is_constant() { return true; }

void TraitConstant::initialize() {
    // do nothing
}

void TraitConstant::generatePhenotypes() {
    // We don't need a GP map for this
    assign(get_dims(), pop.size(), Xinit);
}

/*
void TraitConstant::add_to_checkpoint(Checkpoint &cp) {
    // do nothing
}

void TraitConstant::read_to_checkpoint(Checkpoint &cp, iSimfile &isf) {
    // niente
}

int TraitConstant::resume_from_checkpoint(Checkpoint &cp, int dataIndex){
    return dataIndex;
}
*/

/* class Trait_sample {
    protected:
    int dims;
    int size;
    std::vector<traitType> trait_values;
    public:
    Trait_sample(Trait& T);
    int get_dims() {return dims;}
    int get_size() {return size;}
    void write_to_file(oSimfile &osf);
};*/

Trait_sample::Trait_sample(Trait & T) : Matrix<traitType>(T){
    //assert(dims*size == trait_values.size());
}

//void Trait_sample::write_to_file(oSimfile &osf) {
//    osf.write<size_type>(dims);
//    osf.write<size_type>(size);
//    osf.writeArray<traitType>(&trait_values.at(0), dims*size);
//}

Trait_sample::Trait_sample(iSimfile& isf) : Matrix<traitType>(isf){
//    dims = isf.read<size_type>();
//    size = isf.read<size_type>();
//    isf.readArray<traitType>(trait_values, dims*size);
}
