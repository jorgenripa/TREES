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


#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <functional>

#include "population.h"
#include "mating.h"
#include "sample.h"
#include "randgen.h"
#include "fastmath.h"

Population::Population(ParameterFile& pf) {
    std::string modType;
    std::string modName;

    // Read general parameters:
    gene_tracking = pf.getBool("gene_tracking");
    gene_sampling = pf.getBool("gene_sampling");

    F = pf.getPositiveInt("f");
    
    n0 = (int)pf.getPositiveInt("n_0");
    n = 0;
    age = 0;
    alive.reserve(10000); // Boolean vector, keeping track of who's alive
    alive.assign(0, true); // empty vector
    fitness.reserve(10000);
    fitness.assign(0,1);
    
    // Read genetics model and parameters:
    readGenetics(pf);
    
    // Read traits:
    readTraits(pf, modType, modName);
    
    // read Space module, if specified:
    if (modType=="space") {
        addSpace(modName,pf);
        pf.getStringPair( modType, modName );
    } else {
        space = new NullSpace(*this);
    }
    
    // Read fitness modules
    fitnessList.clear();
    while (modType=="fitness") {
        addFitness(modName,pf);
        pf.getStringPair( modType, modName );
    }
    
    // Mating
    readMating(pf, modType, modName);
    
    
    if (modType.length()>0) {
        std::cout << "Unexpected module at end : " << modType << "\n";
        exit(0);
    }
    
    // Set efficiency flags:
    withinPatchMating = false;
    if (space->getDims()==0) {
        withinPatchMating = true;
    }
    if (space->isDiscrete() && theMatingType==Local && mate_s_space==0.0) {
        withinPatchMating = true;
    }
}

void Population::readGenetics( ParameterFile& pf) {
    std::string geneticsModel = pf.getStringLower("genetics");
    if ( geneticsModel == "continuous_alleles") {
        genetics = new ContinuousAlleles(*this, pf);
    } else if ( geneticsModel == "diallelic" ) {
        genetics = new Diallelic(*this, pf);
    } else {
        std::cout << "Unknown Genetics model!\n";
        exit(0);
    }
}

void Population::readTraits( ParameterFile& pf, std::string& modType, std::string& modName) {
    traitList.clear();
    pf.getStringPair(modType, modName );
    while (modType.substr(0,5)=="trait") {
        if (modType=="trait") {
            addTrait(modName, pf);
        } else if (modType =="trait_constant") {
            addTraitConstant(modName, pf);
        } else {
            std::cout << "Unknown Trait type:" << modType << "\n";
            exit(0);
        }
        pf.getStringPair( modType, modName );
        while (modType=="transform") {
            addTransformToLastTrait(modName,pf);
            pf.getStringPair( modType, modName );
        }
    }
}

void Population::readMating(ParameterFile& pf, std::string& modType, std::string& modName) {
    if (modType=="mating_pool") {
        switch (std::tolower(modName[0])) {
            case 's':
                theMatingType = Selfing;
                break;
            case 'g':
                theMatingType = Random;
                break;
            case 'l':
                theMatingType = Local;
                mate_s_space = pf.getDouble("s_space");
                break;
            default:
                std::cout << "Unknown Mating_pool type : " << modName << "\n";
                exit(0);
                break;
        }
    }else {
        std::cout << "Expected Mating_pool, found : " << modType << "\n";
        exit(0);
    }

    mating_trials = pf.getPositiveInt("mating_trials");
    
    // Read preference modules
    matingPreferenceList.clear();
    pf.getStringPair( modType, modName );
    while (modType=="mating_preference") {
        addPreference(modName,pf);
        pf.getStringPair( modType, modName );
    }
}

void Population::allAlive() {
    alive.assign(n, true);
    somebodyDied = false;
}

void Population::kill(int individual) {
    alive.at(individual) = false;
    somebodyDied = true;
}

void Population::addTrait(std::string& traitName, ParameterFile& pf){
    traitList.push_back(new Trait(traitName, *this, pf));
}
void Population::addTraitConstant(std::string& traitName, ParameterFile& pf){
    traitList.push_back(new TraitConstant(traitName, *this, pf));
}


void Population::addTransformToLastTrait(std::string& tname, ParameterFile& pf) {
    if (tname=="linear") {
        traitList.back()->addTransform(new LinearTransform(pf));
    } else if (tname=="abs") {
        traitList.back()->addTransform(new AbsTransform(pf));
    } else if (tname=="logistic") {
        traitList.back()->addTransform(new LogisticTransform(pf));
    } else if (tname=="normal_deviate") {
        traitList.back()->addTransform(new NormalDeviate(pf));
    } else {
        std::cout << "Unknown transform : " << tname << '\n';
        exit(0);
    }
}


Trait* Population::findTrait(std::string& tname) {
    for (Trait* t :traitList) {
        if (t->getName() == tname) {
            return t;
        }
    }
    std::cout << "Unknown trait: " << tname << '\n';
    exit(0);
}


void Population::addFitness(std::string modName, ParameterFile &pf) {
    //Convert to lower case:
    std::transform(modName.begin(), modName.end(), modName.begin(), ::tolower);

    if (modName=="stabilizing_selection") {
        fitnessList.push_back(new StabilizingSelection(*this, pf));
    } else if (modName=="density_dependence") {
        fitnessList.push_back(new DensityDependence(*this, pf));
    } else if (modName == "resource_landscape") {
        fitnessList.push_back(new ResourceLandscape(*this, pf));
    } else if (modName=="discrete_resources") {
        fitnessList.push_back(new DiscreteResources(*this, pf));
    } else if (modName=="spatial_gradient") {
        fitnessList.push_back(new SpatialGradient(*this, pf));
    } else if (modName=="catastrophes") {
        fitnessList.push_back(new Catastrophe(*this, pf));
//    } else if (modName=="FinchModel") {
//        fitnessList.push_back(new FinchModel(*this, pf));
    } else {
        std::cout << "Unknown Fitness model : " << modName << "\n";
        exit(0);
    }
}

void Population::addPreference(std::string modName, ParameterFile &pf) {
    //Convert to lower case:
    std::transform(modName.begin(), modName.end(), modName.begin(), ::tolower);

    if ( modName == "target_selection") {
        matingPreferenceList.push_back(new Preference(*this, pf));
    } else {
        std::cout << "Unknown Mating_preference model : " << modName << "\n";
        exit(0);
    }
}

void Population::addSpace(std::string modName, ParameterFile &pf) {
    //Convert to lower case:
    std::transform(modName.begin(), modName.end(), modName.begin(), ::tolower);
    if ( modName == "none") {
        space = new NullSpace(*this);
    } else if ( modName == "discrete") {
        space = new DiscreteSpaceImp(*this, pf);
    } else if ( modName == "continuous") {
        space = new ContinuousSpace(*this, pf);
    } else {
        std::cout << "Unknown Space model!\n";
        exit(0);
    }
}

void Population::makeNextGeneration() {
    reproduce(); // mate, produce offspring, replace parents with children
    space->disperse(); // may include dispersal mortality
    // This is to avoid testing for mortality during fitness calculations:
    if (somebodyDied) {
        compactData();
    }
    survive(); // adult survival.
    ++age;
    // Every now and then prune the Allele tables:
    if (gene_tracking && age % 50 == 0) {
        genetics->pruneGeneLists();
    }
}

void Population::reproduce() {
    int nnew = n*F;
    genetics->prepareNewGeneration(nnew);
    space->prepareNewGeneration(nnew);
    std::vector<int> dadList = findDads();
    nnew = 0; // restart counting
    for (int mom=0; mom<n; ++mom) {
        int dad = dadList[mom];
        if (dad>=0) { // there is an accepted mate?
            for (int child=0; child<F; ++child) {
                // recombine from both parents
                genetics->addChild(mom,dad);
                // inherit spatial position (of mom, usually):
                space->addChild(mom, dad);
                ++nnew;
            }
        }
    }
    // Establish next generation:
    n = nnew;
    genetics->nextGeneration();
    space->nextGeneration();
    generatePhenotypes();
    alive.assign(n, true);
    somebodyDied = false;
    fitness.assign(n,1.0);
}

std::vector<int> Population::findDads() {
    std::vector<int> dads(n,-1);
    if (theMatingType == Selfing) {
        for (int mom=0; mom<n; ++mom) {
            dads[mom] = mom;
        }
    } else if (withinPatchMating) {
        DiscreteSpace& theSpace = dynamic_cast<DiscreteSpace&>(getSpace());
        std::vector<int> indList;
        for (int patch=0; patch<theSpace.getPatchCount(); ++patch) {
            if (theSpace.popSizeInPatch(patch)>0) {
                int local_n = theSpace.popSizeInPatch(patch);
                // Generate list of patch inhabitants
                indList.clear();
                for (int i=0; i<size(); ++i) {
                    if (theSpace.getLinearPatch(i) == patch) {
                        indList.push_back(i);
                    }
                }
                // Let each Mom choose from the local list
                for (int momi=0; momi<local_n; ++momi) {
                    int mom = indList[momi];
                    int trialCount = 0;
                    while (trialCount<mating_trials) {
                        // First find a possible candidate
                        int picked = (int)(rand1()*local_n);
                        int candidate = indList[picked];
                        // Next, let mom accept or reject according to her preferences
                        double Pchoose = 1;
                        for (Preference* pp : matingPreferenceList) {
                            Pchoose *= pp->getPartnerWeight(mom, candidate);
                        }
                        if (Pchoose==1 || rand1()<Pchoose) {
                            dads[mom] = candidate;// acceptance!
                            break;
                        }
                        ++trialCount; // rejection, try another mate
                    }
                }
            }
        }
    } else if (space->getDims()==0 || theMatingType == Random) {
        // no space:
        for (int mom=0; mom<n; ++mom) {
            int trialCount = 0;
            while (trialCount<mating_trials) {
                // Pick a random candidate
                int candidate = (int)(rand1()*n);
                // Next, let mom accept or reject according to her preferences
                double Pchoose = 1;
                for (Preference* pp : matingPreferenceList) {
                    Pchoose *= pp->getPartnerWeight(mom, candidate);
                }
                if (Pchoose==1 || rand1()<Pchoose) {
                    dads[mom] = candidate; // acceptance!
                    break; // while loop
                }
                ++trialCount; // rejection, try another mate
            }
        }
    } else {
        // General case
        Fastexp fexp;
        double twoss2 = 2 * mate_s_space * mate_s_space;
        if (/* DISABLES CODE */ (false)) {
            for (int mom=0; mom<n; ++mom) {
                std::vector<double> space_weights;
                space_weights.assign(n, 1);
                // Calculate spatial filter weights:
                for (int male=0; male<n; ++male) {
                    double dist2 = getSpace().getDist2(mom, male);
                    space_weights[male] *= fexp.exp(-dist2/twoss2);
                }
                // calculate cumulative sum:
                std::partial_sum(space_weights.begin(), space_weights.end(), space_weights.begin());
                int trialCount = 0;
                while (trialCount<mating_trials) {
                    // First find a possible candidate
                    // It's a random pick, weighted by spatial filters
                    int candidate = weightedChoiceCumSum(space_weights);
                    // Next, let mom accept or reject according to her preferences
                    double Pchoose = 1;
                    for (Preference* pp : matingPreferenceList) {
                        Pchoose *= pp->getPartnerWeight(mom, candidate);
                    }
                    if (Pchoose==1 || rand1()<Pchoose) {
                        dads[mom] = candidate; // acceptance!
                        break; // while loop
                    }
                    ++trialCount; // rejection, try another mate
                }
            }
        }else { // other solution, considerably faster
            for (int mom=0; mom<n; ++mom) {
                int trialCount = 0;
                while (trialCount<mating_trials) {
                    // First find a possible candidate
                    // It's a random pick, weighted by spatial filters
                    int candidate = -1;
                    while (candidate<0) {
                        int male = int(rand1()*n);
                        double dist2 = getSpace().getDist2(mom, male);
                        if(rand1() < fexp.exp(-dist2/twoss2)) {
                            candidate = male;
                        }
                        
                    }
                    // Next, let mom accept or reject according to her preferences
                    double Pchoose = 1;
                    for (Preference* pp : matingPreferenceList) {
                        Pchoose *= pp->getPartnerWeight(mom, candidate);
                    }
                    if (Pchoose==1 || rand1()<Pchoose) {
                        dads[mom] = candidate; // acceptance!
                        break; // while loop
                    }
                    ++trialCount; // rejection, try another mate
                }
            }
            
        }
    }
    return dads;
}


void Population::generatePhenotypes() {
    for (Trait* tr : traitList) {
        tr->generatePhenotypes();
    }
    fitness.assign(n, 1);
}

void Population::survive() {
    for (Fitness* sp : fitnessList) {
        sp->aggregateFitness(fitness);
    }
    for (int i=0; i<n; ++i) {
        alive.at(i) = alive.at(i) && rand1()<fitness.at(i)/F;
    }
    compactData();
}

void Population::compactData(){
    // compact all data:
    genetics->compactData(alive);
    for( Trait* tr : traitList) {
        tr->compactData(alive);
    }
    space->compactData(alive);
    // compact fitness:
    int iw = 0;
    for (int ir=0; ir<alive.size(); ++ir) {
        if (alive[ir]) {
            fitness[iw]=fitness[ir];
            ++iw;
        }
    }
    n = iw;
    alive.assign(n, true);
    somebodyDied = false;
}

Population::~Population() {
    //std::cout << "killing population\n";
    for( Trait* tr : traitList) {
        delete tr;
    }
    for( Fitness* sp : fitnessList) {
        delete sp;
    }
    delete space;
    for( Preference* pp : matingPreferenceList) {
        delete pp;
    }
    delete genetics;
}

void Population::initialize() {
    // initialize genomes:
    genetics->initialize(n0);
    // set genes to initial values:
    for( Trait* tr : traitList) {
        tr->initialize(n0);
    }
    if (gene_tracking) {
        // create gene tables and assign initial values:
        genetics->initializeTracking(n0);
    }
    // assign initial position:
    space->initialize(n0);
    n = n0;
    alive.assign(n, true);
    generatePhenotypes();
    age = 0;
}

Sample* Population::makeSample() {
    Sample* theSample = new Sample(age,n);
    if (gene_sampling) {
        genetics->addToSample(*theSample);
    }
    for(Trait* tr : traitList) {
        if (!tr->isConstant()) {
            tr->addToSample(*theSample);
        }
    }
    space->addToSample(*theSample);
    if (gene_tracking) {
        genetics->addGeneIDsToSample(*theSample);
    }
    return theSample;
}

