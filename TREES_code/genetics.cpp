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
#include "genotype_phenotype_map.h"
#include "randgen.h"
#include "simfile.h"


Genetics::Genetics( Population& p, ParameterFile& pf) :
pop(p) {
    loci = 0;
    Pmut = pf.getPositiveDouble("p_mutation");

    geneLists.clear();
    if (geneTracking()) {
        G1id.reserve(1,10000);
        G2id.reserve(1,10000);
        newG1id.reserve(1,10000);
        newG2id.reserve(1,10000);
    }
}

Genetics::~Genetics() {}

//bool Genetics::isClonal() { return false; } // default reply

bool Genetics::geneTracking()
     { return pop.isTrackingGenes(); }

void Genetics::initializeTracking() {
    id_type firstId=0;
    for (int l=0; l<loci; ++l) {
        geneLists.push_back(GeneList());
        firstId = geneLists.back().addGene(0,0,getEffect1(0,l));//parent, time, effect
        geneLists.back().getGene(firstId).sampled = true;
    }
    G1id.assign(loci,pop.size(),firstId);
    G2id.assign(loci,pop.size(),firstId);
}

void Genetics::prepareNewGeneration(int size) {
    if (geneTracking()) {
        if ( newG1id.capacity() < loci*size ){
            newG1id.reserve(loci,2*size);
            newG2id.reserve(loci,2*size);
        }
        newG1id.resize(loci,size);
        newG2id.resize(loci,size);
    }
    newCount = 0;
}

void Genetics::nextGeneration() {
    if (geneTracking()) {
        newG1id.resize(loci,newCount); // some matings may fail
        newG2id.resize(loci,newCount); // some matings may fail
        G1id.swap(newG1id);
        G2id.swap(newG2id);
        //newG1id.clear();
        //newG2id.clear();
    }
}

void Genetics::mark_all_living_genes_sampled() {
    // find all current alleles in population,
    // mark them and their ancestors as sampled:
    for (int li=0; li<loci; ++li) {
        GeneList& list = geneLists[li];
        for (int i=0; i<pop.size(); ++i) {
            id_type id = getGene1id(i, li);
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

void Genetics::writeGeneLists(oSimfile &sf) {
    //One list per locus:
    for (GeneList& list : geneLists) {
        list.writeGenes(sf);
    }
}

void Genetics::readGeneLists(iSimfile &sf) {
    //One list per locus
    geneLists.clear();
    geneLists.reserve(loci);
    for (int li=0; li<loci; ++li) {
        geneLists.push_back(GeneList(sf));
    }
}

//void Genetics::addToCheckpoint(Checkpoint &cp) {
//    // save genomes;
//    addGenomesToSample(cp);
//    if (pop.isTrackingGenes()) {
//        addGeneIDsToSample(cp);
//        cp.geneListsCopy = geneLists; // deep copy?
//    }
//}

void Genetics::resumeFromCheckpoint(Genetics_sample* gsample,std::vector<GeneList>& glists) {
    resumeGenomesFromCheckpoint(gsample);
    if (pop.isTrackingGenes()) {
        G1id = gsample->G1id;
        G2id = gsample->G2id;
        geneLists = glists; // deep copy?
    }
}

//void Genetics::readGeneListsToCheckpoint(Checkpoint& cp, iSimfile& isf) {
//    std::vector<GeneList>& gls = cp.geneListsCopy;
//    gls.clear();
//    gls.reserve(loci);
//    for (int li=0; li<loci; ++li) {
//        gls.push_back(GeneList(isf));
//    }
//}
//
void Genetics::pruneGeneLists() {
    // remove all alleles that are not present in the population or have no present descendents
    for(int li=0; li<loci; ++li) {
        GeneList& list = geneLists[li];
        std::vector<bool> alive;
        alive.assign(list.size(), false);
        // find all current alleles in population:
        for (int i=0; i<pop.size(); ++i) {
            id_type id;
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
                id_type parent = list[ai].parent;
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
        G1id.compact_data(alive);
        G2id.compact_data(alive);
    }
}

/*class Genetics_sample {
    protected:
    bool gene_tracking;
    int loci;
    int pop_size;
    std::vector<id_type> G1id;
    std::vector<id_type> G2id;
    Genetics_sample(Genetics& G);
    public:
    virtual void write_to_file(oSimfile& osf);
};*/

Genetics_sample::Genetics_sample() {
    loci = 0;
    pop_size = 0;
    gene_tracking = false;
    G1id.clear();
    G2id.clear();
}

Genetics_sample::Genetics_sample(Genetics& G) {
    loci = G.getLoci();
    pop_size = G.pop.size();
    gene_tracking = G.pop.isTrackingGenes();
    if (gene_tracking) {
        G1id = G.G1id;
        G2id = G.G2id;
        G.mark_all_living_genes_sampled();
    }
}

void Genetics_sample::write_to_file(oSimfile& osf) {
    osf.write<size_type>(loci);
    osf.write<size_type>(pop_size);
    osf.write<char>((char)gene_tracking);
    if (gene_tracking) {
        G1id.write_to_file(osf);
        G2id.write_to_file(osf);
    }
}

Genetics_sample::Genetics_sample(iSimfile& isf) {
    loci = isf.read<size_type>();
    pop_size = isf.read<size_type>();
    gene_tracking = (bool)isf.read<char>();
    if (gene_tracking) {
        G1id.read_from_file(isf);
        G2id.read_from_file(isf);
    }
}

Genetics_sample::~Genetics_sample() {
}

Genetics_sample* Genetics_sample::create_from_file(iSimfile& isf) { // polymorphic read
    Genetics_types type = (Genetics_types)isf.read<char>();
    switch (type) {
        case Genetics_types::Continuous_type: {
            ContinuousAlleles_sample* CAsam = new ContinuousAlleles_sample(isf);
            return CAsam;
            break;
        }
        case Genetics_types::Diallelic_type: {
            Diallelic_sample* Dsam = new Diallelic_sample(isf);
            return Dsam;
            break;
        }
        default:
        std::cout << "Wrong genetics type\n";
        exit(1);
        return NULL;
    }
}

//////////////////////////////////////////////////////////////////////////
// ContinuousAlleles
//////////////////////////////////////////////////////////////////////////

ContinuousAlleles::ContinuousAlleles(Population& pop, ParameterFile& pf) :
Genetics(pop, pf) {
    G1.reserve(1,10000);
    G2.reserve(1,10000);
    newG1.reserve(1,10000);
    newG2.reserve(1,10000);
}

ContinuousAlleles::~ContinuousAlleles() {
}

double ContinuousAlleles::getEffect1(int individual, int locus) {
    return getGene1(individual, locus);
}
double ContinuousAlleles::getEffect2(int individual, int locus) {
    return getGene2(individual, locus);
}

Genotype_Phenotype_map* ContinuousAlleles::create_GP_map(Trait& T, int loci_per_dim, ParameterFile& pf) {
    return new ContinuousAlleles_GP_map(T, *this, loci_per_dim, pf);
}

void ContinuousAlleles::initialize() {
    double P_No_Mutations = pow(1-Pmut,loci);
    P_At_least_one_mut = 1 - P_No_Mutations;

    G1.assign(loci,pop.size(),0);
    G2.assign(loci,pop.size(),0);
    newG1.clear();
    newG2.clear();
}

void ContinuousAlleles::prepareNewGeneration(int size) {
    Genetics::prepareNewGeneration(size);
    if ( newG1.capacity() < loci*size ){
        newG1.reserve(loci,2*size);
        newG2.reserve(loci,2*size);
    }
    newG1.resize(loci,size);
    newG2.resize(loci,size);
}

void ContinuousAlleles::nextGeneration() {
    Genetics::nextGeneration();
    newG1.resize(loci,newCount); // some matings may fail
    newG2.resize(loci,newCount); // some matings may fail
    G1.swap(newG1);
    G2.swap(newG2);
    //newG1.clear();
    //newG2.clear();
}

void ContinuousAlleles::addChild(int mom, int dad) {
    produceGamete(mom,1);
    produceGamete(dad,2);
    ++newCount;
}

void ContinuousAlleles::produceGamete(int parent, int targetHaplo) {
    if (geneTracking()) {
        float* g1 = &getGene1(parent,0);
        float* g2 = &getGene2(parent,0);
        id_type* id1 = &getGene1id(parent,0);
        id_type* id2 = &getGene2id(parent,0);
        float* target;
        id_type* idTarget;
        if (targetHaplo==1) {
            target = &newG1(0,newCount);
            idTarget = &newG1id(0,newCount);
        } else {
            target = &newG2(0,newCount);
            idTarget = &newG2id(0,newCount);
        }
        bitGenerator bG;
        float* target_it = target;
        id_type* idTarget_it = idTarget;
        for (int a=0; a<loci; ++a) {
            if (bG.nextBit()) { //} haplo ==1) {
                *target_it = g1[a];
                *idTarget_it = id1[a]; // parallel copying of allele id
            } else {
                *target_it = g2[a];
                *idTarget_it = id2[a]; // parallel copying of allele id
            }
            ++target_it;
            ++idTarget_it;
        }
        // Possibly mutate:
        while (rand1() < P_At_least_one_mut) {
            int mutlocus = rand1()*loci;
            target[mutlocus] += doubleExp(1); // mutations have variance 1
            idTarget[mutlocus] = geneLists[mutlocus].addGene(idTarget[mutlocus], pop.getAge(), target[mutlocus]);
        }
    } else { // no gene tracking
        float* g1 = &getGene1(parent,0);
        float* g2 = &getGene2(parent,0);
        float* target;
        if (targetHaplo==1) {
            target = &newG1(0,newCount);
        } else {
            target = &newG2(0,newCount);
        }
        
        bitGenerator bG;
        float* target_it = target;
        for (int a=0; a<loci; ++a) {
            if (bG.nextBit()) { //(haplo ==1) {
                *target_it = g1[a];
            } else {
                *target_it = g2[a];
            }
            ++target_it;
        }
        // Possibly mutate:
        while (rand1() < P_At_least_one_mut) {
            int mutlocus = rand1()*loci;
            target[mutlocus] += doubleExp(1); // mutations have variance 1
        }
        
    }
}

void ContinuousAlleles::compactData(std::vector<bool>& alive) {
    Genetics::compactData(alive);
    G1.compact_data(alive);
    G2.compact_data(alive);
}

void ContinuousAlleles::resumeGenomesFromCheckpoint(Genetics_sample* gsample) {
    initialize();
    ContinuousAlleles_sample * CA_sample = dynamic_cast<ContinuousAlleles_sample*>(gsample);
    G1 = CA_sample->G1;
    G2 = CA_sample->G2;
}

int ContinuousAlleles::add_loci(int L) {
    int start_locus = loci;
    loci += L;
    return start_locus;
}

Genetics_sample* ContinuousAlleles::make_sample() {
    return new ContinuousAlleles_sample(*this);
}

/*class ContinuousAlleles_sample : public Genetics_sample {
    std::vector<ContinuousAlleles::geneType> G1;
    std::vector<ContinuousAlleles::geneType> G2;
    
    public:
    ContinuousAlleles_sample(ContinuousAlleles& G);
    virtual void write_to_file(oSimfile& osf);
};*/

ContinuousAlleles_sample::ContinuousAlleles_sample(ContinuousAlleles& G) : Genetics_sample(G) {
    G1 = G.G1;
    G2 = G.G2;
}

void ContinuousAlleles_sample::write_to_file(oSimfile& osf) {
    // First specify type, convenient for reading:
    osf.write<char>((char)Genetics_types::Continuous_type);
    // Let base class write basic data and gene tracking data, if applicable
    Genetics_sample::write_to_file(osf);
    G1.write_to_file(osf);
    G2.write_to_file(osf);
}

ContinuousAlleles_sample::ContinuousAlleles_sample(iSimfile& isf) :
Genetics_sample(isf) {
    G1.read_from_file(isf);
    G2.read_from_file(isf);
}

// This virtual destructor is essential.
ContinuousAlleles_sample::~ContinuousAlleles_sample() {
}

///////////////////////////////////////////////////////////////////
//  Diallelic
///////////////////////////////////////////////////////////////////

Diallelic::Diallelic(Population& pop, ParameterFile& pf) :
Genetics(pop, pf) {
    G1.reserve(1,10000);
    G2.reserve(1,10000);
    newG1.reserve(1,10000);
    newG2.reserve(1,10000);
}

Diallelic::~Diallelic() {
}

double Diallelic::getEffect1(int individual, int locus) {
    return double(getGene1(individual, locus))/2.0;
}
double Diallelic::getEffect2(int individual, int locus) {
    return double(getGene2(individual, locus))/2.0;
}

Genotype_Phenotype_map* Diallelic::create_GP_map(Trait& T, int loci_per_dim, ParameterFile& pf) {
    return new Diallelic_GP_map(T, *this, loci_per_dim, pf);
}

void Diallelic::initialize() {
    double P_No_Mutations = pow(1-Pmut,loci);
    P_At_least_one_mut = 1 - P_No_Mutations;
    
    G1.assign(loci,pop.size(),-1);
    G2.assign(loci,pop.size(),-1);
    newG1.clear();
    newG2.clear();
}

void Diallelic::prepareNewGeneration(int size) {
    Genetics::prepareNewGeneration(size);
    if ( newG1.capacity() < loci*size ){
        newG1.reserve(loci,2*size);
        newG2.reserve(loci,2*size);
        //G1.reserve(2*loci*size);
        //G2.reserve(2*loci*size);
    }
    newG1.resize(loci,size);
    newG2.resize(loci,size);
    //void* c2 = (void*)&newG1.getX().at(0);
}

void Diallelic::nextGeneration() {
    Genetics::nextGeneration(); // deal with gene tracking
    newG1.resize(loci,newCount); // some matings may fail
    newG2.resize(loci,newCount); // some matings may fail
    G1.swap(newG1);
    G2.swap(newG2);
    //newG1.clear();
    //newG2.clear();
}

void Diallelic::addChild(int mom, int dad) {
    produceGamete(mom,1);
    produceGamete(dad,2);
    ++newCount;
}

void Diallelic::produceGamete(int parent, int targetHaplo) {
    if (geneTracking()) {
        geneType* g1 = &getGene1(parent,0);
        geneType* g2 = &getGene2(parent,0);
        id_type* id1 = &getGene1id(parent,0);
        id_type* id2 = &getGene2id(parent,0);
        geneType* target;
        id_type* idTarget;
        if (targetHaplo==1) {
            target = &newG1(0,newCount); //.at(loci*newCount);
            idTarget = &newG1id(0,newCount);
        } else {
            target = &newG2(0,newCount);
            idTarget = &newG2id(0,newCount);
        }
        bitGenerator bG;
        //int haplo = 1 + int(rand1()*2); // copy from parent haplotype 1 or 2
        geneType* target_it = target;
        id_type* idTarget_it = idTarget;
        for (int a=0; a<loci; ++a) {
            if (bG.nextBit()) { // (haplo ==1) {
                *target_it = g1[a];
                *idTarget_it = id1[a]; // parallel copying of allele id
            } else {
                *target_it = g2[a];
                *idTarget_it = id2[a]; // parallel copying of allele id
            }
            ++target_it;
            ++idTarget_it;
            /*if( bG.nextBit() ) { // rand1()<0.5){
                haplo = 3 - haplo; // changes 1 to 2 and vice versa
            }*/
        }
        // Possibly mutate:
        while (rand1() < P_At_least_one_mut) {
            int mutlocus = rand1()*loci;
            target[mutlocus] = -target[mutlocus]; // mutations are always of mean size=1
            double targetEffect = double(target[mutlocus])/2.0;
            idTarget[mutlocus] = geneLists[mutlocus].addGene(idTarget[mutlocus], pop.getAge(), targetEffect);
        }
    } else { // no gene tracking
        geneType* g1 = &getGene1(parent,0);
        geneType* g2 = &getGene2(parent,0);
        geneType* target;
        if (targetHaplo==1) {
            target = &newG1(0,newCount);
        } else {
            target = &newG2(0,newCount);
        }
        
        bitGenerator bG;
        //int haplo = 1 + int(rand1()*2); // copy from parent haplotype 1 or 2
        geneType* target_it = target;
        //std::memcpy(target_it,g2,loci*sizeof(geneType));
        for (int a=0; a<loci; ++a) {
            if (bG.nextBit()) {// (haplo ==1) {
                *target_it = g1[a];
            } else {
                *target_it = g2[a];
            }
            ++target_it;
            /*if( bG.nextBit()){
                haplo = 3 - haplo; // changes 1 to 2 and vice versa
            }*/
        }
        // Possibly mutate:
        while (rand1() < P_At_least_one_mut) {
            int mutlocus = rand1()*loci;
            target[mutlocus] = -target[mutlocus];
        }
        
    }
}

void Diallelic::compactData(std::vector<bool>& alive) {
    Genetics::compactData(alive);
    // compact genome:
    G1.compact_data(alive);
    G2.compact_data(alive);
}

void Diallelic::resumeGenomesFromCheckpoint(Genetics_sample* gsample) {
    Diallelic_sample* D_sample = dynamic_cast<Diallelic_sample*>(gsample);
    
    // copy stored bitvectors to genes
    // (0/1) correspond to (-1/+1)
    initialize(); // sets all genes to -1
    for (int li=0; li<loci; ++li) {
        for (int ind=0; ind<pop.size(); ++ind) {
            if (D_sample->G1.get_bit(li + ind*loci)) {
                G1(li,ind) = +1;
            }
            if (D_sample->G2.get_bit(li + ind*loci)) {
                G2(li,ind) = +1;
            }
        }
    }
}

int Diallelic::add_loci(int L) {
    int start_locus = loci;
    loci += L;
    return start_locus;
}

Genetics_sample* Diallelic::make_sample() {
    return new Diallelic_sample(*this);
}

Diallelic_sample::Diallelic_sample(Diallelic& G) : Genetics_sample(G) {
    G1.reserve(loci*pop_size);
    G2.reserve(loci*pop_size);
    for (int ind=0; ind<pop_size; ++ind) {
        for (int locus=0; locus<loci; ++locus) {
            G1.push_back(G.getGene1(ind, locus)>0);
            G2.push_back(G.getGene2(ind, locus)>0);
        }
    }
}

void Diallelic_sample::write_to_file(oSimfile& osf) {
    osf.write<char>((char)Genetics_types::Diallelic_type);
    // First let base class write basic data and gene tracking data, if applicable
    Genetics_sample::write_to_file(osf);
    G1.write_to_file(osf);
    G2.write_to_file(osf);
}

Diallelic_sample::Diallelic_sample(iSimfile& isf) :
Genetics_sample(isf) {
    G1.read_from_file(isf);
    G2.read_from_file(isf);
}

// This virtual destructor is essential.
Diallelic_sample::~Diallelic_sample() {
}

/////////////////////////////////////////
// Omnigenic genetics model
/////////////////////////////////////////
Omnigenic::Omnigenic(Population& pop, ParameterFile& pf) :
ContinuousAlleles(pop, pf) {
    loci = pf.getPositiveInt("loci"); // This is only set once for this genetics module
}

Genotype_Phenotype_map* Omnigenic::create_GP_map(Trait& T, int loci_per_dim, ParameterFile& pf) {
    return new Omnigenic_GP_map(T, *this, loci_per_dim, pf);
}
    
