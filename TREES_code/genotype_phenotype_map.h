//
//  genotype_phenotype_map.hpp
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

#ifndef genotype_phenotype_map_h
#define genotype_phenotype_map_h

#include <vector>
#include "types.h"
#include "parameterFile.h"

class Genetics;
class Trait;
class Transform;
class Diallelic;
class ContinuousAlleles;
class Omnigenic;
class Checkpoint;
class iSimfile;
class oSimfile;

class Genotype_Phenotype_map {
protected:
    Trait& T;
    Genotype_Phenotype_map(Trait& theTrait);
    std::vector<Transform*> transforms;
    void read_all_transforms(ParameterFile& pf);
    double do_all_transforms(double y);
public:
    virtual ~Genotype_Phenotype_map();
    virtual void generate_phenotypes()=0;
    virtual void initialize(traitType X_init)=0;
    virtual std::vector<double>* make_checkpoint_data();
    virtual void resume_from_checkpoint_data(std::vector<double>* data);
    int get_dims();
    int pop_size();
    Trait& get_trait();
};

class Diallelic_GP_map : public Genotype_Phenotype_map {
    Diallelic& G;
    int start_locus;
    int loci_per_dim;
public:
    Diallelic_GP_map(Trait& theT, Diallelic& theG, int loci_per_dim, ParameterFile& pf);
    virtual ~Diallelic_GP_map();
    virtual void generate_phenotypes();
    virtual void initialize(traitType X_init);
    int get_loci_per_dim() { return loci_per_dim;}
};

class ContinuousAlleles_GP_map : public Genotype_Phenotype_map {
    ContinuousAlleles& G;
    int start_locus;
    int loci_per_dim;
public:
    ContinuousAlleles_GP_map(Trait& theT, ContinuousAlleles& theG, int loci_per_dim, ParameterFile& pf);
    virtual ~ContinuousAlleles_GP_map();
    virtual void generate_phenotypes();
    virtual void initialize(traitType X_init);
};

class Omnigenic_GP_map : public Genotype_Phenotype_map {
    Omnigenic& G;
    std::vector<double> weights;
    double X_init;
public:
    Omnigenic_GP_map(Trait& theT, Omnigenic& theG, int loci_per_dim, ParameterFile& pf);
    virtual void generate_phenotypes();
    virtual void initialize(traitType X_init);
    virtual std::vector<double>* make_checkpoint_data();
    virtual void resume_from_checkpoint_data(std::vector<double>* data);
};

#endif /* genotype_phenotype_map_hpp */
