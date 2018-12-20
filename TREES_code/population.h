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
#include "sample.h"
#include "types.h"

////////////////////////////////
// Population class
// The Population class is where the current state of the population is stored,
// and where all ecological and genetical parameters are kept.
////////////////////////////////
class Population {
protected:
    int n; //population size (sometimes including dead)
    int n0; // initial population size
    int F; // constant fecundity
    int mating_trials; // maximal number of males a female can reject
    timeType age;
    bool gene_tracking;
    bool gene_sampling;
    bool withinPatchMating; // used for efficiency
    Genetics* genetics;
    std::vector<Fitness*> fitnessList;
    MatingType theMatingType;
    double mate_s_space;
    std::vector<Preference*> matingPreferenceList;
    std::vector<int> findDads();
    std::vector<Trait*>  traitList; // a list of traits
    Space* space;
    std::vector<bool>   alive;
    std::vector<double> fitness;
    bool somebodyDied;
    void allAlive();
    void readGenetics( ParameterFile& pf);
    void readTraits( ParameterFile& pf, std::string& modType, std::string& modName);
    void addFitness(std::string modName, ParameterFile& pf);
    void readMating(ParameterFile& pf, std::string& modType, std::string& modName);
    void addPreference(std::string modName, ParameterFile& pf);
    void addSpace(std::string modName, ParameterFile& pf);

    void compactData();
    void generatePhenotypes();
    void addTrait(std::string& traitName, ParameterFile& pf);
    void addTraitConstant(std::string& traitName, ParameterFile& pf);
    void addTransformToLastTrait(std::string& tname, ParameterFile& pf);
    double& getTrait(int individual, int trait) { return traitList.at(trait)->traitValue(individual); }
    std::vector<Trait*>& getTraits() {return traitList;}
public:
    Population(ParameterFile& pf);
    ~Population();
    // Added as a quick fix, remove in Super4:
    void initialize();
    Sample* makeSample();
    void makeNextGeneration(); //
    void reproduce();
    void survive();
    void kill(int individual);
    int size() { return n;}
    Genetics& getGenetics() {return *genetics;}
    Space& getSpace() {return *space;}
    int getF() { return F;}
    Trait* findTrait(std::string& name);
    timeType getAge() { return age;}
    bool isTrackingGenes() { return gene_tracking; }
};

#endif /* defined(__Species__population__) */
