//
//  genotype_phenotype_map.cpp
//  TREES
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
#include "genotype_phenotype_map.h"
#include "population.h"
#include "genetics.h"
#include "trait.h"
#include "transform.h"
#include "randgen.h"

Genotype_Phenotype_map::Genotype_Phenotype_map(Trait& theTrait) : T(theTrait) {
    transforms.clear();
}

Genotype_Phenotype_map::~Genotype_Phenotype_map() {
    for (Transform* t : transforms) {
        delete t;
    }
}

double Genotype_Phenotype_map::do_all_transforms(double y) {
    for (Transform* t : transforms) {
        y = t->transform_single(y);
    }
    return y;
}

int Genotype_Phenotype_map::get_dims() {
    return T.get_dims();
}

int Genotype_Phenotype_map::pop_size() {
    return T.pop.size();
}

Trait& Genotype_Phenotype_map::get_trait() {
    return T;
}

void Genotype_Phenotype_map::read_all_transforms(ParameterFile& pf) {
    while (pf.get_next_name()=="transform") {
        std::string tname = pf.get_next_value_string_lower();
        if (tname=="linear") {
            transforms.push_back(new LinearTransform(pf));
        } else if (tname=="abs") {
            transforms.push_back(new AbsTransform(pf));
        } else if (tname=="logistic") {
            transforms.push_back(new LogisticTransform(pf));
        } else if (tname=="normal_deviate") {
            transforms.push_back(new NormalDeviate(pf));
        } else if (tname=="range") {
            transforms.push_back(new Range(pf,*this));
        } else {
            std::cout << "Unknown transform : " << tname << '\n';
            exit(1);
        }
    }
}

std::vector<double>* Genotype_Phenotype_map::make_checkpoint_data() {
    // Default is to store nothing
    // All necessary information is normally taken from the parameter file.
    return NULL;
}

void Genotype_Phenotype_map::resume_from_checkpoint_data(std::vector<double> *data){
    // default is to do nothing
    // All necessary information is normally taken from the parameter file.
}


Diallelic_GP_map::Diallelic_GP_map(Trait& theT, Diallelic& theG, int loci_per_dim, ParameterFile& pf) : Genotype_Phenotype_map(theT), G(theG) {
    this->loci_per_dim = loci_per_dim;
    start_locus = G.add_loci(loci_per_dim * theT.get_dims());
    read_all_transforms(pf);
}

Diallelic_GP_map::~Diallelic_GP_map() {
}

void Diallelic_GP_map::generate_phenotypes() {
    T.assign(T.get_dims(), pop_size(), 0);
    for (int ind=0; ind<pop_size(); ++ind) {
        for (int dim=0; dim<T.get_dims(); ++dim) {
            double y = 0;
            Diallelic::geneType* g1p = &G.getGene1(ind, start_locus + dim*loci_per_dim);
            Diallelic::geneType* g2p = &G.getGene2(ind, start_locus + dim*loci_per_dim);
            for (int li=0; li<loci_per_dim; ++li) {
                y += *g1p; //G.getGene1(ind, start_locus + dim*loci_per_dim + li);
                y += *g2p; //G.getGene2(ind, start_locus + dim*loci_per_dim + li);
                ++g1p;
                ++g2p;
            }
            // Diallelic genes are coded as +-1, but
            // the effects should actually be +-1/2:
            y = y/2;
            T.traitValue(ind,dim) = do_all_transforms(y);
        }
    }
}

void Diallelic_GP_map::initialize(traitType X_init) {
    int plusCount = (int)((X_init+loci_per_dim)/2.0 + 0.5);
    if (plusCount<0 || plusCount>loci_per_dim) {
        std::cout << "Illegal initial trait value " << X_init << " of trait " << T.getName() << ".\n";
        std::cout << "The allowed range is [-L, L], where L is the number of loci per dimension.\n";
        std::cout << "Note that the initial value is an untransformed value.\n";
        exit(1);
    }
    // set the first loci to +1, the rest remain -1
    for (int ind=0; ind<pop_size(); ++ind) {
        for (int dim=0; dim<T.get_dims(); ++dim) {
            for (int l=0; l<plusCount; ++l) {
                G.getGene1(ind, start_locus + dim*loci_per_dim + l) = 1;
                G.getGene2(ind, start_locus + dim*loci_per_dim + l) = 1;
            }
        }
    }
}

//////////////////////////////////////////////////
// ContinuousAlleles_GP_map
//////////////////////////////////////////////////

ContinuousAlleles_GP_map::ContinuousAlleles_GP_map(Trait& theT, ContinuousAlleles& theG, int loci_per_dim, ParameterFile& pf) : Genotype_Phenotype_map(theT), G(theG) {
    this->loci_per_dim = loci_per_dim;
    start_locus = G.add_loci(loci_per_dim * theT.get_dims());
    read_all_transforms(pf);
}

ContinuousAlleles_GP_map::~ContinuousAlleles_GP_map() {
}

void ContinuousAlleles_GP_map::generate_phenotypes() {
    T.assign(T.get_dims(), pop_size(), 0);
    for (int ind=0; ind<pop_size(); ++ind) {
        for (int dim=0; dim<T.get_dims(); ++dim) {
            double y = 0;
            for (int li=0; li<loci_per_dim; ++li) {
                y += G.getGene1(ind, start_locus + dim*loci_per_dim + li);
                y += G.getGene2(ind, start_locus + dim*loci_per_dim + li);
            }
            T.traitValue(ind,dim) = do_all_transforms(y);
        }
    }
}

void ContinuousAlleles_GP_map::initialize(traitType X_init) {
    // sets first locus to initial value, the rest remain zero
    for (int ind=0; ind<pop_size(); ++ind) {
        for(int d=0; d<T.get_dims(); ++d) {
            G.getGene1(ind, start_locus + d*loci_per_dim) = X_init/2.0;
            G.getGene2(ind, start_locus + d*loci_per_dim) = X_init/2.0;
        }
    }
}


//////////////////////////////////////////////////
// Omnigenic_GP_map
//////////////////////////////////////////////////

Omnigenic_GP_map::Omnigenic_GP_map(Trait& theT, Omnigenic& theG, int loci_per_dim, ParameterFile& pf) : Genotype_Phenotype_map(theT), G(theG) {
    // choose weights:
    int L = G.getLoci();
    int D = T.get_dims();
    weights.assign(L*D, 0);
    for (int di=0; di<D; ++di) {
        double w2sum = 0;
        for (int li=0; li<L; ++li) {
            double w = doubleExp(1);
            weights.at(di*L + li) = w;
            w2sum += w*w;
        }
        // rescale weights such that added mutational variance is proportional to loci_per_dim
        for (int li=0; li<L; ++li) {
            weights.at(di*L + li) *= sqrt(loci_per_dim/w2sum);
        }
    }
    read_all_transforms(pf);
}

void Omnigenic_GP_map::generate_phenotypes() {
    int L = G.getLoci();
    T.assign(T.get_dims(), pop_size(), 0);
    for (int ind=0; ind<pop_size(); ++ind) {
        for (int dim=0; dim<T.get_dims(); ++dim) {
            double y = X_init;
            for (int li=0; li<G.getLoci(); ++li) {
                y += G.getGene1(ind, li)*weights.at(dim*L + li);
                y += G.getGene2(ind, li)*weights.at(dim*L + li);;
            }
            T.traitValue(ind,dim) = do_all_transforms(y);
        }
    }
}

void Omnigenic_GP_map::initialize(traitType X_init) {
    // The X_init is part of the mapping here instead of an initial genetic state
    // See generate_phenotypes()
    this->X_init = X_init;
}

std::vector<double>* Omnigenic_GP_map::make_checkpoint_data() {
    // return a copy of weights:
    return new std::vector<double>(weights);
}

void Omnigenic_GP_map::resume_from_checkpoint_data(std::vector<double>* data) {
    weights = *data;
}

/*
void Omnigenic_GP_map::add_to_checkpoint(Checkpoint &cp) {
    // First add X_init as a vector with a single element:
    std::vector<double> X_init_vector(1, X_init);
    cp.addData(new XData<double>(X_init_vector));
    cp.addData(new XData<double>(weights));
}

void Omnigenic_GP_map::read_to_checkpoint(Checkpoint &cp, iSimfile &isf) {
    cp.addData(new XData<double>(isf,1)); // Read X_init parameter
    cp.addData(new XData<double>(isf,G.getLoci()*T.get_dims())); // Read weights
}

int Omnigenic_GP_map::resume_from_checkpoint(Checkpoint &cp, int dataIndex) {
    X_init = (dynamic_cast<XData<double>&>(cp.getData(dataIndex++))).getData().at(0);
    weights = (dynamic_cast<XData<double>&>(cp.getData(dataIndex++))).getData();
    return dataIndex;
}
*/

