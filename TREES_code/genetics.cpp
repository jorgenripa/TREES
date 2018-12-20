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


#include <cmath>
#include <iostream>
#include <string.h>
#include "stdlib.h"

#include "genetics.h"
#include "population.h"
#include "sample.h"
#include "randgen.h"
#include "simfile.h"


Genetics::Genetics( Population& p, ParameterFile& pf) :
pop(p) {
    loci = 0;
    Pmut = pf.getPositiveDouble("p_mutation");

    geneLists.clear();
    if (geneTracking()) {
        G1id.reserve(10000);
        G2id.reserve(10000);
        newG1id.reserve(10000);
        newG2id.reserve(10000);
    }
}

Genetics::~Genetics() {}

//bool Genetics::isClonal() { return false; } // default reply

bool Genetics::geneTracking()
     { return pop.isTrackingGenes(); }

void Genetics::initializeTracking(int n0) {
    idType firstId;
    for (int l=0; l<loci; ++l) {
        geneLists.push_back(GeneList());
        firstId = geneLists.back().addGene(0,0,getEffect1(0,l));//parent, time, effect
        geneLists.back().getGene(firstId).sampled = true;
    }
    G1id.assign(n0*loci,firstId);
    G2id.assign(n0*loci,firstId);
}

void Genetics::prepareNewGeneration(int size) {
    if (geneTracking()) {
        newG1id.assign(loci*size,0);
        newG2id.assign(loci*size,0);
    }
    newCount = 0;
}

void Genetics::nextGeneration() {
    if (geneTracking()) {
        newG1id.resize(newCount*loci); // some matings may fail
        newG2id.resize(newCount*loci); // some matings may fail
        G1id = newG1id;
        G2id = newG2id;
        newG1id.clear();
        newG2id.clear();
    }
}

void Genetics::addGeneIDsToSample(Sample& s) {
    s.addData(new IdData(G1id));
    s.addData(new IdData(G2id));
    // find all current alleles in population,
    // mark them and their ancestors as sampled:
    for (int li=0; li<loci; ++li) {
        GeneList& list = geneLists[li];
        for (int i=0; i<pop.size(); ++i) {
            idType id = getGene1id(i, li);
            Gene* gp = &list.getGene(id);
            while (!gp->sampled) {
                gp->sampled = true;
                gp = &list.getGene(gp->parent);
            }
            id = getGene2id(i, li);
            gp = &list.getGene(id);
            while (!gp->sampled) {
                gp->sampled = true;
                gp = &list.getGene(gp->parent);
            }
        }
    }
}

void Genetics::pruneGeneLists() {
    // remove all alleles that are not present in the population, have no present descendents, and have never been sampled
    for(int li=0; li<loci; ++li) {
        GeneList& list = geneLists[li];
        std::vector<bool> alive;
        alive.assign(list.size(), false);
        // find all current alleles in population:
        for (int i=0; i<pop.size(); ++i) {
            idType id;
            id = getGene1id(i, li);
            alive.at(list.getGeneIndex(id)) = true;
            id = getGene2id(i, li);
            alive.at(list.getGeneIndex(id)) = true;
        }
        // assign deathtime to recently extinct alleles:
        for (int ai=0; ai<list.size(); ++ai) {
            if (!alive[ai] && list[ai].deathTime==-1) {
                list[ai].deathTime = pop.getAge();
            }
        }
        // make sure their ancestors are not removed:
        for (int ai=0; ai<list.size(); ++ai) {
            if (alive[ai]) {
                idType parent = list[ai].parent;
                while (parent>0) {
                    int parent_i = list.getGeneIndex(parent);
                    if (!alive.at(parent_i)) {
                        alive[parent_i] = true;
                        parent = list[parent_i].parent;
                    } else {
                        parent = 0; // no need to track further back
                    }
                }
            }
        }
        //  prune all dead and non-sampled genes:
        list.pruneGenes(alive);
    }
}

void Genetics::compactData(std::vector<bool>& alive) {
    if (geneTracking()) {
        // compact allele id information:
        int iwrite=0; // position for write
        int iread; // position for read
        for (iread=0; iread<pop.size(); ++iread) {
            if (alive.at(iread)) {
                if(iread>iwrite) {
                    for (int li=0; li<loci; ++li) {
                        getGene1id(iwrite,li) = getGene1id(iread,li);
                        getGene2id(iwrite,li) = getGene2id(iread,li);
                    }
                }
                ++iwrite;
            }
        }
        G1id.resize(iwrite*loci);
        G2id.resize(iwrite*loci);
    }
}

void Genetics::saveGeneLists(Simfile &sf) {
    //One list per locus:
    for (GeneList& list : geneLists) {
        sf.write((int)list.size());
        for (Gene& g : list) {
            sf << g;
        }
    }
}

ContinuousAlleles::ContinuousAlleles(Population& pop, ParameterFile& pf) :
Genetics(pop, pf) {
    G1.reserve(10000);
    G2.reserve(10000);
    newG1.reserve(10000);
    newG2.reserve(10000);
}

ContinuousAlleles::~ContinuousAlleles() {
}

double ContinuousAlleles::getEffect1(int individual, int locus) {
    return getGene1(individual, locus);
}
double ContinuousAlleles::getEffect2(int individual, int locus) {
    return getGene2(individual, locus);
}

void ContinuousAlleles::initialize(int n0) {
    double P_No_Mutations = pow(1-Pmut,loci);
    P_At_least_one_mut = 1 - P_No_Mutations;

    G1.assign(loci*n0,0);
    G2.assign(loci*n0,0);
    newG1.clear();
    newG2.clear();
}

void ContinuousAlleles::setInitialValue(double value, int n, int startLocus, int dims, int lociPerDim ) {
    // sets first locus to initial value, the rest remain zero
    for (int ind=0; ind<n; ++ind) {
        for(int d=0; d<dims; ++d) {
            getGene1(ind, startLocus + d*lociPerDim) = value/2.0;
            getGene2(ind, startLocus + d*lociPerDim) = value/2.0;
        }
    }
}


void ContinuousAlleles::prepareNewGeneration(int size) {
    Genetics::prepareNewGeneration(size);
    newG1.assign(loci*size,0);
    newG2.assign(loci*size,0);
}

void ContinuousAlleles::nextGeneration() {
    Genetics::nextGeneration();
    newG1.resize(newCount*loci); // some matings may fail
    newG2.resize(newCount*loci); // some matings may fail
    G1 = newG1;
    G2 = newG2;
    newG1.clear();
    newG2.clear();
}

void ContinuousAlleles::addChild(int mom, int dad) {
    produceGamete(mom,1);
    produceGamete(dad,2);
    ++newCount;
}

void ContinuousAlleles::produceGamete(int parent, int targetHaplo) {
    if (geneTracking()) {
        double* g1 = &getGene1(parent,0);
        double* g2 = &getGene2(parent,0);
        idType* id1 = &getGene1id(parent,0);
        idType* id2 = &getGene2id(parent,0);
        double* target;
        idType* idTarget;
        if (targetHaplo==1) {
            target = &newG1.at(loci*newCount);
            idTarget = &newG1id.at(loci*newCount);
        } else {
            target = &newG2.at(loci*newCount);
            idTarget = &newG2id.at(loci*newCount);
        }
        bitGenerator bG;
        int haplo = 1 + int(rand1()*2); // copy from parent haplotype 1 or 2
        double* target_it = target;
        idType* idTarget_it = idTarget;
        for (int a=0; a<loci; ++a) {
            if (haplo ==1) {
                *target_it = g1[a];
                *idTarget_it = id1[a]; // parallel copying of allele id
            } else {
                *target_it = g2[a];
                *idTarget_it = id2[a]; // parallel copying of allele id
            }
            ++target_it;
            ++idTarget_it;
            if( bG.nextBit() ) { // rand1()<0.5){
                haplo = 3 - haplo; // changes 1 to 2 and vice versa
            }
            
        }
        // Possibly mutate:
        while (rand1() < P_At_least_one_mut) {
            int mutlocus = rand1()*loci;
            target[mutlocus] += doubleExp(1); // mutations are always of mean size=1
            idTarget[mutlocus] = geneLists[mutlocus].addGene(idTarget[mutlocus], pop.getAge(), target[mutlocus]);
            //assert(nG[mutlocus]>-1 && nG[mutlocus]<1);
        }
    } else { // no gene tracking
        double* g1 = &getGene1(parent,0);
        double* g2 = &getGene2(parent,0);
        double* target;
        if (targetHaplo==1) {
            target = &newG1.at(loci*newCount);
        } else {
            target = &newG2.at(loci*newCount);
        }
        
        bitGenerator bG;
        int haplo = 1 + int(rand1()*2); // copy from parent haplotype 1 or 2
        double* target_it = target;
        for (int a=0; a<loci; ++a) {
            if (haplo ==1) {
                *target_it = g1[a];
            } else {
                *target_it = g2[a];
            }
            ++target_it;
            //assert(nG[a]>-1 && nG[a]<1);
            if( bG.nextBit()){
                haplo = 3 - haplo; // changes 1 to 2 and vice versa
            }
        }
        // Possibly mutate:
        while (rand1() < P_At_least_one_mut) {
            int mutlocus = rand1()*loci;
            target[mutlocus] += randn(); // mutations are always of SD=1
            //assert(nG[mutlocus]>-1 && nG[mutlocus]<1);
        }
        
    }
}

double ContinuousAlleles::getGeneSum(int individual, int startlocus, int endlocus) {
    // no gene interactions:
    double sum = 0;
    double* g1 = &getGene1(individual,startlocus);
    double* g2 = &getGene2(individual,startlocus);
    for (int li=startlocus; li<endlocus; ++li) {
        sum += *g1 + *g2;
        ++g1;
        ++g2;
    }
    return sum;
    //return getGene1(individual, locus) + getGene2(individual, locus);
}

void ContinuousAlleles::compactData(std::vector<bool>& alive) {
    Genetics::compactData(alive);
    // compact genome:
    int iwrite=0; // position for write
    int iread; // position for read
    for (iread=0; iread<pop.size(); ++iread) {
        if (alive.at(iread)) {
            if(iread>iwrite) {
                memcpy(&getGene1(iwrite, 0), &getGene1(iread, 0), loci*sizeof(G1[0]));
                memcpy(&getGene2(iwrite, 0), &getGene2(iread, 0), loci*sizeof(G1[0]));
/*                for (int a=0; a<loci; ++a) {
                    getGene1(iwrite,a) = getGene1(iread,a);
                    getGene2(iwrite,a) = getGene2(iread,a);
                }
 */
            }
            ++iwrite;
        }
    }
    G1.resize(iwrite*loci);
    G2.resize(iwrite*loci);
}

void ContinuousAlleles::addToSample(Sample& s) {
    s.addData(new FloatData(G1));
    s.addData(new FloatData(G2));
}


////////////////////////////////

Diallelic::Diallelic(Population& pop, ParameterFile& pf) :
Genetics(pop, pf) {
    G1.reserve(10000);
    G2.reserve(10000);
    newG1.reserve(10000);
    newG2.reserve(10000);
}

Diallelic::~Diallelic() {
}

double Diallelic::getEffect1(int individual, int locus) {
    return (double)getGene1(individual, locus);
}
double Diallelic::getEffect2(int individual, int locus) {
    return (double)getGene2(individual, locus);
}


void Diallelic::initialize(int n0) {
    double P_No_Mutations = pow(1-Pmut,loci);
    P_At_least_one_mut = 1 - P_No_Mutations;
    
    G1.assign(loci*n0,-1);
    G2.assign(loci*n0,-1);
    newG1.clear();
    newG2.clear();

}

void Diallelic::setInitialValue(double value, int n, int startLocus, int dims, int lociPerDim ) {
    //int traitLoci = dims*lociPerDim;
    int plusCount = (int)((value+lociPerDim)/2.0);
    if (plusCount<0 || plusCount>lociPerDim) {
        std::cout << "Illegal initial trait value : " << value << '\n';
        exit(0);
    }
    // set the first loci to +1, the rest remain -1
    for (int ind=0; ind<n; ++ind) {
        for (int d=0; d<dims; ++d) {
            for (int l=0; l<plusCount; ++l) {
                getGene1(ind, startLocus + d*lociPerDim + l) = 1;
                getGene2(ind, startLocus + d*lociPerDim + l) = 1;
            }
        }
    }
}

void Diallelic::prepareNewGeneration(int size) {
    Genetics::prepareNewGeneration(size);
    newG1.assign(loci*size,0);
    newG2.assign(loci*size,0);
}

void Diallelic::nextGeneration() {
    Genetics::nextGeneration(); // deal with gene tracking
    newG1.resize(newCount*loci); // some matings may fail
    newG2.resize(newCount*loci); // some matings may fail
    G1 = newG1;
    G2 = newG2;
    newG1.clear();
    newG2.clear();
}

void Diallelic::addChild(int mom, int dad) {
    produceGamete(mom,1);
    produceGamete(dad,2);
    ++newCount;
}

void Diallelic::produceGamete(int parent, int targetHaplo) {
    if (geneTracking()) {
        int* g1 = &getGene1(parent,0);
        int* g2 = &getGene2(parent,0);
        idType* id1 = &getGene1id(parent,0);
        idType* id2 = &getGene2id(parent,0);
        int* target;
        idType* idTarget;
        if (targetHaplo==1) {
            target = &newG1.at(loci*newCount);
            idTarget = &newG1id.at(loci*newCount);
        } else {
            target = &newG2.at(loci*newCount);
            idTarget = &newG2id.at(loci*newCount);
        }
        bitGenerator bG;
        int haplo = 1 + int(rand1()*2); // copy from parent haplotype 1 or 2
        int* target_it = target;
        idType* idTarget_it = idTarget;
        for (int a=0; a<loci; ++a) {
            if (haplo ==1) {
                *target_it = g1[a];
                *idTarget_it = id1[a]; // parallel copying of allele id
            } else {
                *target_it = g2[a];
                *idTarget_it = id2[a]; // parallel copying of allele id
            }
            ++target_it;
            ++idTarget_it;
            if( bG.nextBit() ) { // rand1()<0.5){
                haplo = 3 - haplo; // changes 1 to 2 and vice versa
            }
        }
        // Possibly mutate:
        while (rand1() < P_At_least_one_mut) {
            int mutlocus = rand1()*loci;
            target[mutlocus] = -target[mutlocus]; // mutations are always of mean size=1
            idTarget[mutlocus] = geneLists[mutlocus].addGene(idTarget[mutlocus], pop.getAge(), target[mutlocus]);
        }
    } else { // no gene tracking
        int* g1 = &getGene1(parent,0);
        int* g2 = &getGene2(parent,0);
        int* target;
        if (targetHaplo==1) {
            target = &newG1.at(loci*newCount);
        } else {
            target = &newG2.at(loci*newCount);
        }
        
        bitGenerator bG;
        int haplo = 1 + int(rand1()*2); // copy from parent haplotype 1 or 2
        int* target_it = target;
        for (int a=0; a<loci; ++a) {
            if (haplo ==1) {
                *target_it = g1[a];
            } else {
                *target_it = g2[a];
            }
            ++target_it;
            if( bG.nextBit()){
                haplo = 3 - haplo; // changes 1 to 2 and vice versa
            }
        }
        // Possibly mutate:
        while (rand1() < P_At_least_one_mut) {
            int mutlocus = rand1()*loci;
            target[mutlocus] = -target[mutlocus];
        }
        
    }
}

double Diallelic::getGeneSum(int individual, int startlocus, int endlocus) {
    // no gene interactions:
    double sum = 0;
    int* g1 = &getGene1(individual,startlocus);
    int* g2 = &getGene2(individual,startlocus);
    for (int li=startlocus; li<endlocus; ++li) {
        sum += *g1 + *g2;
        ++g1;
        ++g2;
    }
    return sum/2; // Divide by 2 here, such that the effects are +/- 1/2
}

void Diallelic::compactData(std::vector<bool>& alive) {
    Genetics::compactData(alive);
    // compact genome:
    int iwrite=0; // position for write
    int iread; // position for read
    for (iread=0; iread<pop.size(); ++iread) {
        if (alive.at(iread)) {
            if(iread>iwrite) {
                memcpy(&getGene1(iwrite, 0), &getGene1(iread, 0), loci*sizeof(G1[0]));
                memcpy(&getGene2(iwrite, 0), &getGene2(iread, 0), loci*sizeof(G1[0]));
            }
            ++iwrite;
        }
    }
    G1.resize(iwrite*loci);
    G2.resize(iwrite*loci);
}

void Diallelic::addToSample(Sample& s) {
    s.addData(new IntData(G1));
    s.addData(new IntData(G2));
}
