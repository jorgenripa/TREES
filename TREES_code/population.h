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


#ifndef __Species__population__
#define __Species__population__

#include <stdio.h>
#include <vector>

#include "genetics.h"
#include "fitness.h"
#include "mating.h"
#include "space.h"
#include "trait.h"
#include "parameterFile.h"
#include "types.h"

class Population_checkpoint;

////////////////////////////////
// Population class
// The Population class is where the current state of the population is stored,
// and where all ecological and genetical parameters are kept.
////////////////////////////////
class Population {
    friend class Population_sample;
    friend class Population_checkpoint;
    friend class Microsample;
protected:
// Constant simulation parameters:
    int F; // constant fecundity
    int mating_trials; // maximal number of males a female can reject
    bool gene_tracking;
    bool gene_sampling;
    bool withinPatchMating; // used for efficiency
    MatingType theMatingType;
    double mate_s_space;
    int n0; // initial population size
    
// Model structures:
    Genetics* genetics;
    Space* space;
    std::vector<Preference*> matingPreferenceList;
    std::vector<Trait*>  traitList; // a list of traits
    std::vector<Fitness*> fitnessList;

// Dynamic parameters:
    int n; //population size (sometimes including dead)
    time_type age;
    std::vector<bool>   alive;
    std::vector<double> fitness;
    bool somebodyDied;
    
    std::vector<int> findDads();
    void allAlive();
    void readGenetics( ParameterFile& pf);
    void readTraits( ParameterFile& pf);
    void addFitness(ParameterFile& pf);
    void readMating(ParameterFile& pf);
    void addPreference(ParameterFile& pf);
    void addSpace(ParameterFile& pf);

    void compactData();
    void generatePhenotypes();
    void addTrait(std::string traitName, ParameterFile& pf);
    void addTraitConstant(std::string traitName, ParameterFile& pf);
    //void addTransformToLastTrait(std::string tname, ParameterFile& pf);
    traitType& getTrait(int individual, int trait) { return traitList.at(trait)->traitValue(individual); }
    std::vector<Trait*>& getTraits() {return traitList;}
    int calc_total_data_dimensions();


public:
    Population(ParameterFile& pf);
    ~Population();
    void initialize();
    void resumeAtCheckpoint(Population_checkpoint& cp);
    // Running:
    void makeNextGeneration(); //
    void reproduce();
    void survive();
    void kill(int individual);
    // Info:
    int size() { return n;}
    Genetics& getGenetics() {return *genetics;}
    Space& getSpace() {return *space;}
    int getF() { return F;}
    Trait* findTrait(std::string& name);
    time_type getAge() { return age;}
    bool isTrackingGenes() { return gene_tracking; }
};

class Sample_base {
protected:
    time_type generation;
    double cputime; // in seconds
    int pop_size;
    
    void set(time_type gen, double cpu, int n)
        { generation=gen; cputime=cpu; pop_size=n;}
    Sample_base(time_type generation, double cputime, int n)
        { set(generation,cputime,n);}
    Sample_base()
        {set(0,0,0);}
    
public:
    time_type get_generation() { return generation; }
    double get_cputime() { return cputime;}
    int get_pop_size() {return pop_size;}
    virtual ~Sample_base();
};


class Population_sample : public Sample_base {
    protected:
    Genetics_sample* genes_sample;
    std::vector<Trait_sample> traits;
    Space_sample* space_sample;
    public:
    Population_sample(Population& pop, double cputime);
    Population_sample(iSimfile& isf);
    virtual ~Population_sample();
    void write_to_file(oSimfile& osf);
};

class Population_checkpoint : public Sample_base {
    friend class Population;
    // stores genes, genelists and space, and possible trait parameters
    // (trait values can be reconstructed from genes)
    seed_type seed;
    Genetics_sample* genetics;
    std::vector<GeneList> geneListsCopy;
    Space_sample* space;
    // A vector of vectors (one per non-constant trait, or none):
    std::vector<std::vector<double>*> GP_map_data; // possible trait parameters (see Omnigenic model)
    
    public:
    Population_checkpoint(Population& pop, seed_type cpseed, double cputime);
    Population_checkpoint(iSimfile& isf);
    ~Population_checkpoint();
    void write_to_file(oSimfile& osf);
    seed_type get_seed() { return seed;}
    std::vector<double>* get_map_data(int mi) { return GP_map_data.at(mi);}
};

class Microsample : public Sample_base {
    // optionally stores means, variances and covariances of all traits and spatial positions
    protected:
    void calc_means(Population& pop);
    void calc_variances(Population& pop);
    void calc_covariances(Population& pop);
    
    public:
    char option;
    std::vector<traitType> means;
    std::vector<traitType> covariances; // either just variances or all covariances
    Microsample(Population & pop, char option, double cputime);
    Microsample(iSimfile& isf);
    void write_to_file(oSimfile& osf);
};

#endif /* defined(__Species__population__) */
