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


#ifndef __Species__trait__
#define __Species__trait__

#include <stdio.h>
#include <vector>
#include <string>

#include "transform.h"
#include "parameterFile.h"
#include "types.h" // defines traitType
#include "matrix.h"

class Population;
class Genetics;
class Genotype_Phenotype_map;
class Checkpoint;
class iSimfile;

/////////////////////////////
// Trait class
// Contains all trait values
// Has a genotype-phenotype map
/////////////////////////////
class Trait : public Matrix<traitType> {
    friend class Population;
    friend class Genotype_Phenotype_map;
    friend class Diallelic_GP_map;
    friend class ContinuousAlleles_GP_map;
    friend class Omnigenic_GP_map;
    friend class Trait_sample;
protected:
    Trait(std::string& name, Population& p);
    std::string name;
    Population& pop;
    Genetics& genes;
    Genotype_Phenotype_map* GP_map;
    void set_dims(int d) { M = d;}
    traitType Xinit; // initial value
public:
    Trait(std::string& name, Population& p, ParameterFile& pf);
    ~Trait();
    virtual bool is_constant();
    virtual void initialize();
    virtual void generatePhenotypes();
    void compactData(std::vector<bool>& alive) { compact_data(alive);}
    traitType& traitValue(int individual, int dim=0);
    std::string& getName() { return name; }
    int get_dims() { return M; }
    Genotype_Phenotype_map* get_GP_map() { return GP_map;}
//    virtual void add_to_checkpoint(Checkpoint& cp);
//    virtual void read_to_checkpoint(Checkpoint& cp, iSimfile& isf);
//    virtual int resume_from_checkpoint(Checkpoint& cp, int dataIndex);
};

class TraitConstant : public Trait {
public:
    TraitConstant(std::string& name, Population& p, ParameterFile& pf);
    virtual bool is_constant();
    virtual void initialize();
    virtual void generatePhenotypes();
//    virtual void add_to_checkpoint(Checkpoint& cp);
//    virtual void read_to_checkpoint(Checkpoint& cp, iSimfile& isf);
//    virtual int resume_from_checkpoint(Checkpoint& cp, int dataIndex);
};

class Trait_sample : public Matrix<traitType>{
protected:
public:
    Trait_sample(Trait& T);
    Trait_sample(iSimfile& isf);
    int get_dims() {return M;}
    int get_size() {return N;}
    //void write_to_file(oSimfile &osf);
};

#endif /* defined(__Species__trait__) */
