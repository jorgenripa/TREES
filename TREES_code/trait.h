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

class Population;
class Genetics;
class Sample;

/////////////////////////////
// Trait class
// A trait is really a genotype-phenotype map
// The Trait class is friends with the Genetics class
/////////////////////////////
class Trait {
protected:
    Trait(std::string& name, Population& p);
    std::string name;
    Population& pop;
    Genetics& genes;
    int dims; // number of dimensions
    std::vector<double> X; // phenotypic data of entire population
    double Xinit; // initial value
    typedef std::vector<double>::iterator Xiter;
    int startLocus;
    int lociPerDim; // loci per dimension
    std::vector<Transform*> transforms;
public:
    Trait(std::string& name, Population& p, ParameterFile& pf);
    ~Trait();
    void addTransform( Transform* t);
    virtual void initialize(int n0);
    virtual void generatePhenotypes();
    void compactData(std::vector<bool>& alive);
    double& traitValue(int individual, int dim=0);
    //    {return X[individual];} // not worth it!
    //std::vector<double>& getX() { return X;}
    std::string& getName() { return name; }
    int getDims() { return dims; }
    void addToSample(Sample& s);
    virtual bool isConstant();
};

class TraitConstant : public Trait {
public:
    TraitConstant(std::string& name, Population& p, ParameterFile& pf);
    virtual void initialize(int n0);
    virtual void generatePhenotypes();
    virtual bool isConstant();
};

/*
class MultiTraitIterator {
protected:
    int traitCount; // height of virtual array
    double** traitp; // list of pointers to trait values
public:
    MultiTraitIterator(std::vector<Trait>& traitList, int startpos=0);
    MultiTraitIterator(std::vector<Trait*>& traitList, int startpos=0);
    ~MultiTraitIterator();
    inline double& operator [](int t) {return *traitp[t];}
    inline MultiTraitIterator& operator ++() {
        for (int t=0; t<traitCount; ++t) {
            ++traitp[t];
        }
        return *this;
    }
    
    //bool atEnd();
};
*/

#endif /* defined(__Species__trait__) */
