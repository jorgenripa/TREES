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


#ifndef __Species__genetics__
#define __Species__genetics__

#include <stdio.h>
#include "parameterFile.h"
#include "sample.h"
#include "geneTracking.h"
#include "types.h"

class Population;

///////////////////////////////////////
// Genetics
// Base class for genetics modules
// Handles gene tacking, if activated
////////////////////////////////////////
class Genetics {
    friend class Trait;
    friend class TraitConstant;
protected:
    Population& pop;
    int loci; // the number of loci
    double Pmut;
    int newCount; // keeps count of next generation
    virtual double getEffect1(int individual, int locus)=0;
    virtual double getEffect2(int individual, int locus)=0;

    // Tracking:
    std::vector<GeneList> geneLists; // Gene tables for tracking, one per locus
    void addGene(int locus, idType parent, double effect, timeType time);
    std::vector<idType> G1id, G2id; // mirrors of G1 and G2
    std::vector<idType> newG1id, newG2id; // mirrors of newG1 and newG2
    inline idType& getGene1id(int individual, int locus) {
        return G1id.at(individual*loci + locus);
    }
    inline idType& getGene2id(int individual, int locus) {
        return G2id.at(individual*loci + locus);
    }

    // protected constructor:
    Genetics( Population& p, ParameterFile& pf);
public:
    virtual ~Genetics();
    //virtual bool isClonal(); // returns false by default
    void initializeTracking(int n0);
    virtual void initialize(int n0)=0;
    virtual void setInitialValue(double value, int n, int startLocus, int dims, int lociPerDim )=0;
    virtual void prepareNewGeneration(int size);
    virtual void nextGeneration();
    virtual void addChild(int mom, int dad)=0;
    virtual double getGeneSum(int individual, int startlocus, int endlocus)=0; // endlocus is exclusive!
    virtual void compactData(std::vector<bool>& alive);
    virtual void addToSample(Sample& s)=0;
    virtual void checkGenes()=0;

    void addGeneIDsToSample(Sample& s);
    int getLoci() { return loci;}
    void pruneGeneLists();
    bool geneTracking();
    void saveGeneLists(Simfile& sf);
};

// Continuum of alleles genetics:
class ContinuousAlleles : public Genetics {
    friend class Trait;
public:
    ContinuousAlleles(Population& pop, ParameterFile& pf);
    virtual ~ContinuousAlleles();
    virtual void initialize(int n0);
    virtual void setInitialValue(double value, int n, int startLocus, int dims, int lociPerDim );
    virtual void prepareNewGeneration(int size);
    virtual void nextGeneration();
    virtual void addChild(int mom, int dad);
    virtual double getGeneSum(int individual, int startlocus, int endlocus); // endlocus is exclusive!
    virtual void compactData(std::vector<bool>& alive);
    virtual void addToSample(Sample& s);
    virtual void checkGenes();

protected:
    virtual double getEffect1(int individual, int locus);
    virtual double getEffect2(int individual, int locus);
    double P_At_least_one_mut;
    std::vector<double> G1, G2;
    //int Gsize;
    std::vector<double> newG1, newG2;
    //int newSize;
    inline double& getGene1(int individual, int locus) {
        return G1.at(individual*loci + locus);
    }
    inline double& getGene2(int individual, int locus) {
        return G2.at(individual*loci + locus);
    }
    void produceGamete(int parent, int targetHaplo);
};


// Diallelic genetics (0/1 alleles):
class Diallelic : public Genetics {
    friend class Trait;
public:
    Diallelic(Population& pop, ParameterFile& pf);
    virtual ~Diallelic();
    virtual void initialize(int n0);
    virtual void setInitialValue(double value, int n, int startLocus, int dims, int lociPerDim );
    virtual void prepareNewGeneration(int size);
    virtual void nextGeneration();
    virtual void addChild(int mom, int dad);
    virtual double getGeneSum(int individual, int startlocus, int endlocus); // endlocus is exclusive!
    virtual void compactData(std::vector<bool>& alive);
    virtual void addToSample(Sample& s);
    
protected:
    virtual double getEffect1(int individual, int locus);
    virtual double getEffect2(int individual, int locus);
    double P_At_least_one_mut;
    std::vector<int> G1, G2;
    std::vector<int> newG1, newG2;
    inline int& getGene1(int individual, int locus) {
        return G1.at(individual*loci + locus);
    }
    inline int& getGene2(int individual, int locus) {
        return G2.at(individual*loci + locus);
    }
    void produceGamete(int parent, int targetHaplo);
    // Debugging:
    virtual void checkGenes();
    void checkChild(int ci, int haplo);
    void checkAllChildren();
};

#endif /* defined(__Species__genetics__) */
