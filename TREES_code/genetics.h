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
#include "geneTracking.h"
#include "types.h"
#include "matrix.h"
#include "bit_vector.h"

class Population;
class Genotype_Phenotype_map;
class Trait;
class ParameterFile;

///////////////////////////////////////
// Genetics
// Base class for genetics modules
// Handles gene tacking, if activated
////////////////////////////////////////
class Genetics_sample;

class Genetics {
//    friend class Trait;
//    friend class TraitConstant;
    friend class Genetics_sample;
protected:
    Population& pop;
    int loci; // the number of loci
    double Pmut;
    int newCount; // keeps count of next generation
    virtual double getEffect1(int individual, int locus)=0;
    virtual double getEffect2(int individual, int locus)=0;

    // Tracking:
    std::vector<GeneList> geneLists; // Gene tables for tracking, one per locus
    void addGene(int locus, id_type parent, double effect, time_type time);
    Matrix<id_type> G1id, G2id; // mirrors of G1 and G2
    Matrix<id_type> newG1id, newG2id; // mirrors of newG1 and newG2
    inline id_type& getGene1id(int individual, int locus) {
        return G1id(locus,individual);
    }
    inline id_type& getGene2id(int individual, int locus) {
        return G2id(locus,individual);
    }

    // protected constructor:
    Genetics( Population& p, ParameterFile& pf);
public:
    virtual ~Genetics();
    //virtual bool isClonal(); // returns false by default
    virtual Genotype_Phenotype_map* create_GP_map(Trait& T, int loci_per_dim, ParameterFile& pf)=0;
    void initializeTracking();
    virtual void initialize()=0;
    //virtual void setInitialValue(double value, int n, int startLocus, int dims, int lociPerDim )=0;
    virtual void prepareNewGeneration(int size);
    virtual void nextGeneration();
    virtual void addChild(int mom, int dad)=0;
    //virtual double getGeneSum(int individual, int startlocus, int endlocus)=0; // endlocus is exclusive!
    virtual void compactData(std::vector<bool>& alive);
    virtual void resumeGenomesFromCheckpoint(Genetics_sample* gsample)=0;
    virtual Genetics_sample* make_sample()=0;
    void resumeFromCheckpoint(Genetics_sample* gsample, std::vector<GeneList>& glists);

    int getLoci() { return loci;}
    void pruneGeneLists();
    bool geneTracking();
    void writeGeneLists(oSimfile& sf);
    void readGeneLists(iSimfile& sf);
    void mark_all_living_genes_sampled();
    const std::vector<GeneList>& getGeneLists() { return geneLists; }
};

// Enumeration used for polymorphic file read/write
enum Genetics_types {
    Continuous_type,
    Diallelic_type
};

class Genetics_sample {
    friend class Genetics;
    protected:
    bool gene_tracking;
    int loci;
    int pop_size;
    Matrix<id_type> G1id;
    Matrix<id_type> G2id;
    Genetics_sample();
    Genetics_sample(Genetics& G);
    Genetics_sample(iSimfile& isf); // base class read

    public:
    virtual ~Genetics_sample();
    virtual void write_to_file(oSimfile& osf);
    static Genetics_sample* create_from_file(iSimfile& isf); // polymorphic read
};

// Continuum of alleles genetics:
class ContinuousAlleles : public Genetics {
//    friend class Trait;
    friend class ContinuousAlleles_GP_map;
    friend class ContinuousAlleles_sample;
protected:
    // The local coding is float, to save space
    typedef float geneType;
public:
    ContinuousAlleles(Population& pop, ParameterFile& pf);
    virtual ~ContinuousAlleles();
    virtual Genotype_Phenotype_map* create_GP_map(Trait& T, int loci_per_dim, ParameterFile& pf);
    virtual void initialize();
    //virtual void setInitialValue(double value, int n, int startLocus, int dims, int lociPerDim );
    virtual void prepareNewGeneration(int size);
    virtual void nextGeneration();
    virtual void addChild(int mom, int dad);
    //virtual double getGeneSum(int individual, int startlocus, int endlocus); // endlocus is exclusive!
    virtual void compactData(std::vector<bool>& alive);
    virtual void resumeGenomesFromCheckpoint(Genetics_sample* gsample);
    virtual Genetics_sample* make_sample();

    int add_loci(int L);
protected:
    virtual double getEffect1(int individual, int locus);
    virtual double getEffect2(int individual, int locus);
    double P_At_least_one_mut;
    Matrix<geneType> G1, G2;
    //int Gsize;
    Matrix<geneType> newG1, newG2;
    //int newSize;
    inline geneType& getGene1(int individual, int locus) {
        return G1(locus,individual);
    }
    inline geneType& getGene2(int individual, int locus) {
        return G2(locus,individual);
    }
    void produceGamete(int parent, int targetHaplo);
};

class ContinuousAlleles_sample : public Genetics_sample {
    friend class ContinuousAlleles;
    protected:
    Matrix<ContinuousAlleles::geneType> G1;
    Matrix<ContinuousAlleles::geneType> G2;

    public:
    ContinuousAlleles_sample();
    ContinuousAlleles_sample(ContinuousAlleles& G);
    ContinuousAlleles_sample(iSimfile& isf);
    virtual ~ContinuousAlleles_sample();
    virtual void write_to_file(oSimfile& osf);
};

// Diallelic genetics (+-0.5 alleles):
class Diallelic : public Genetics {
    friend class Diallelic_GP_map;
    friend class Diallelic_sample;
protected:
    // The local coding is +-1, although the effects are actually +-0.5.
    // This is adjusted in 'getEffect1/2' and 'getGeneSum'.
    typedef signed char geneType;
public:
    Diallelic(Population& pop, ParameterFile& pf);
    virtual ~Diallelic();
    virtual Genotype_Phenotype_map* create_GP_map(Trait& T, int loci_per_dim, ParameterFile& pf);
    virtual void initialize();
    //virtual void setInitialValue(double value, int n, int startLocus, int dims, int lociPerDim );
    virtual void prepareNewGeneration(int size);
    virtual void nextGeneration();
    virtual void addChild(int mom, int dad);
    //virtual double getGeneSum(int individual, int startlocus, int endlocus); // endlocus is exclusive!
    virtual void compactData(std::vector<bool>& alive);
    virtual void resumeGenomesFromCheckpoint(Genetics_sample* gsample);
    virtual Genetics_sample* make_sample();
    int add_loci(int L);
protected:
    virtual double getEffect1(int individual, int locus);
    virtual double getEffect2(int individual, int locus);
    double P_At_least_one_mut;
    
    Matrix<geneType> G1, G2;
    Matrix<geneType> newG1, newG2;
    inline geneType& getGene1(int individual, int locus) {
        return G1(locus,individual);
    }
    inline geneType& getGene2(int individual, int locus) {
        return G2(locus,individual);
    }
    void produceGamete(int parent, int targetHaplo);
};

class Diallelic_sample : public Genetics_sample {
    friend class Diallelic;
    bit_vector G1; // compressed format of 0/1 bits
    bit_vector G2;
    
    public:
    Diallelic_sample();
    Diallelic_sample(Diallelic& G);
    Diallelic_sample(iSimfile& isf);
    virtual ~Diallelic_sample();
    virtual void write_to_file(oSimfile& osf);
};

// Omnigenic model, complete pleitropy:
class Omnigenic : public ContinuousAlleles {
    //    friend class Trait;
    friend class Omnigenic_GP_map;
public:
    Omnigenic(Population& pop, ParameterFile& pf);
    virtual Genotype_Phenotype_map* create_GP_map(Trait& T, int loci_per_dim, ParameterFile& pf);
};

#endif /* defined(__Species__genetics__) */
