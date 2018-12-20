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


#include "trait.h"
#include "population.h"
#include "sample.h"
#include "randgen.h"

// protected constructor, used by subclasses:
Trait::Trait(std::string& tname, Population& p) :
name(tname), pop(p), genes(p.getGenetics()) {
    dims = 0;
    lociPerDim = 0;
    Xinit = 0;
    startLocus = 0;
}

Trait::Trait(std::string& tname, Population& p, ParameterFile& pf) :
name(tname), pop(p), genes(p.getGenetics()) {
    
    // Make genomic space for this trait:
    dims = pf.getPositiveInt("dimensions");
    lociPerDim = pf.getPositiveInt("loci_per_dim");
    Xinit =pf.getDouble("initial_value");
    startLocus = genes.loci;
    genes.loci += lociPerDim*dims;
    
    // Reserve memory space for phenotypes (this will be expanded when necessary):
    X.reserve(10000*dims);
    // clear transform list:
    transforms.clear();
}

bool Trait::isConstant() { return false; }

void Trait::initialize(int n0) {
    genes.setInitialValue(Xinit,n0,startLocus,dims,lociPerDim);
}

Trait::~Trait() {
    for (int ti=0; ti<transforms.size(); ++ti) {
        delete transforms.at(ti);
    }
}

void Trait::addTransform(Transform *t) {
    transforms.push_back(t);
}

void Trait::generatePhenotypes() {
    X.assign(pop.size()*dims, 0);
    if (lociPerDim>0) {
        for (int i=0; i<pop.size(); ++i) {
            for (int d=0; d<dims; ++d) {
                traitValue(i,d) = genes.getGeneSum(i, startLocus + d*lociPerDim, startLocus + (d+1)*lociPerDim);
            }
        }
    }
    for (int ti=0; ti<transforms.size(); ++ti) {
        transforms[ti]->transform(&X[0], X.size());
    }
}

double& Trait::traitValue(int individual, int dim) {
    return X.at(individual*dims + dim);
}

void Trait::compactData(std::vector<bool>& alive) {
    // compact arrays:
    int iw = 0; // write index
    for (int ir=0; ir<pop.size(); ++ir) { // ir = read index
        if (alive[ir]) {
            if (ir>iw) {
                for (int d=0; d<dims; ++d) {
                    traitValue(iw,d) = traitValue(ir,d);
                }
            }
            ++iw;
        }
    }
    X.resize(iw*dims);
}

void Trait::addToSample(Sample& s) {
    s.addData(new FloatData(X));
}

TraitConstant::TraitConstant(std::string& name, Population& p, ParameterFile& pf):
Trait(name,p) {
    
    // A constant trait is implemented as a trait with zero loci:
    dims = pf.getPositiveInt("dimensions");
    lociPerDim = 0;
    Xinit =pf.getDouble("initial_value");
    startLocus = genes.loci;
    
    // Reserve memory space for phenotypes (this will be expanded when necessary):
    X.reserve(10000*dims);
    // clear transform list:
    transforms.clear();
}

bool TraitConstant::isConstant() { return true; }

void TraitConstant::initialize(int n0) {
    // do nothing
}

void TraitConstant::generatePhenotypes() {
    X.assign(pop.size()*dims, Xinit);
}

