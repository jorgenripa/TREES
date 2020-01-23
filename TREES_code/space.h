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


#ifndef __Species__space__
#define __Species__space__

#include <stdio.h>
#include <vector>
#include <cmath>

#include "parameterFile.h"
#include "trait.h"
#include "matrix.h"

/////////////////////////////
// Space modules
// A Space object holds the spatial as well as dispersal parameters.
// I keeps track of the positions of all individuals and manages dispersal.
/////////////////////////////

class Population;
class Space_sample;

// The space interface:
class Space {
    friend class Space_sample;
protected:
    Population& pop;
    int Dims;
public:
    Space(Population& p);
    virtual ~Space();
    int getDims() { return Dims; }
    virtual void initialize()=0;
    virtual void disperse()=0;
    virtual void prepareNewGeneration(int size)=0;
    virtual void nextGeneration()=0;
    virtual void addChild(int mom, int dad)=0;
    virtual void compactData(std::vector<bool>& alive)=0;
    virtual Space_sample* make_sample()=0;
    virtual void resumeFromCheckpoint(Space_sample* S_sample)=0;
    virtual double getPosition(int individual, int dimension)=0;
    virtual double getDist2(int ind1, int ind2)=0; // squared distance
    virtual double get_mean(int dim)=0;
    virtual double get_variance(int dim, double mean)=0;
    virtual bool isDiscrete()=0;
};

enum Space_types {
    Discrete_space_type,
    Continuous_space_type
};

class Space_sample {
protected:
    int pop_size;
    int dims;
    Space_sample(Space& S);
    Space_sample(iSimfile& isf);
public:
    virtual ~Space_sample();
    virtual void write_to_file(oSimfile& osf);
    static Space_sample* create_from_file(iSimfile& isf);
};

// The DiscreteSpace interface
class Discrete_space : public Space {
    friend class Discrete_space_sample;
public:
    Discrete_space(Population& p);
    //virtual ~DiscreteSpace()=0;
    // new functions:
    virtual int getLinearPatch(int individual)=0;
    virtual int popSizeInPatch(int linearPatch)=0;
    virtual int getPatchCount()=0;
    // From Space interface:
    virtual void initialize()=0;
    virtual void disperse()=0;
    virtual void prepareNewGeneration(int size)=0;
    virtual void nextGeneration()=0;
    virtual void addChild(int mom, int dad)=0;
    virtual void compactData(std::vector<bool>& alive)=0;
    virtual void resumeFromCheckpoint(Space_sample* S_sample)=0;
    virtual double getPosition(int individual, int dim)=0;
    virtual double getDist2(int ind1, int ind2)=0; // squared distance
    virtual double get_mean(int dim)=0;
    virtual double get_variance(int dim, double mean)=0;
    virtual bool isDiscrete();
    virtual Space_sample* make_sample();
};

class Discrete_space_sample : public Space_sample {
    friend class Discrete_space_imp;
protected:
    std::vector<int> linear_patches;
public:
    Discrete_space_sample(Discrete_space& S);
    Discrete_space_sample(iSimfile& isf);
    virtual ~Discrete_space_sample();
    virtual void write_to_file(oSimfile& osf);
};

// A discrete space implementation
class Discrete_space_imp : public Discrete_space {
protected:
    Matrix<int> patches;
    Matrix<int> newPatches;
    std::vector<int> patchPopSizes;
    std::vector<int> linearPatches; // one per individual
    time_type generationLastCount;
    int length; // size of space = length^Dims
    int initialPatch; // starting patch for all individuals (in all dimensions)
    Trait* PDisp;
    enum DispersalType {Global, Neighbour, Distance};
    DispersalType dispType;
    double dispdist; // only used if dispType is Distance
    double Pjump; // only used if dispType is Distance
    enum Boundary {Circular, Reflective, Absorbing};
    Boundary boundaryCondition;
    int treatBoundaries(int pos, int individual);
    // inline for speed, may be negative:
    inline int getDistance(int pos1, int pos2) {
        int dx = pos1-pos2;
        if (boundaryCondition==Circular) {
            if (dx>length/2) { // (integer division)
                return length-dx;
            } else if (dx<-length/2) {
                return length+dx;
            } else
                return dx;
        } else {
            return dx;
        }
    }
    int calcLinearPatch(int individual);
    std::vector<int> calc_coords_from_linear(int linear);
    void assignLinearPatches();
public:
    Discrete_space_imp(Population& p, ParameterFile& pf);
    virtual ~Discrete_space_imp();
    int& getPatch(int individual, int dim) {return patches(dim,individual);}
    // DiscreteSpace interface:
    virtual int getLinearPatch(int individual);
    virtual int popSizeInPatch(int linearPatch);
    virtual int getPatchCount() { return std::pow(length,Dims); }
    // Space interface:
    virtual void initialize();
    virtual void disperse();
    virtual void prepareNewGeneration(int size);
    virtual void nextGeneration();
    virtual void addChild(int mom, int dad);
    virtual void compactData(std::vector<bool>& alive);
    virtual void resumeFromCheckpoint(Space_sample* S_sample);
    virtual double getPosition(int individual, int dim);
    virtual double getDist2(int ind1, int ind2); // squared distance
    virtual double get_mean(int dim);
    virtual double get_variance(int dim, double mean);
};

// The NullSpace class is a DiscreteSpace implementation
//   with no stored positions.
// All individuals are regarded as positioned in patch 0.
class Null_space : public Discrete_space {
public:
    Null_space(Population& p);
    virtual ~Null_space();

    virtual int getLinearPatch(int individual);
    virtual int popSizeInPatch(int linearPatch);
    virtual int getPatchCount();
    // From Space interface:
    virtual void initialize();
    virtual void disperse();
    virtual void prepareNewGeneration(int size);
    virtual void nextGeneration();
    virtual void addChild(int mom, int dad);
    virtual void compactData(std::vector<bool>& alive);
    virtual void resumeFromCheckpoint(Space_sample* S_s);
    virtual double getPosition(int individual, int dimension);
    virtual double getDist2(int ind1, int ind2); // squared distance
    virtual double get_mean(int dim);
    virtual double get_variance(int dim, double mean);
};


class Continuous_space : public Space {
    friend class Continuous_space_sample;
protected:
    typedef float positionType;
    Matrix<positionType> pos;
    Matrix<positionType> newPos;
    double initialPosition;  // starting position for all individuals (in all dimensions)
    Trait* PDisp;
    positionType dispdist; // mean dispersal distance
    positionType maxPos; // size in each dimension
    enum Boundary {Circular, Reflective, Absorbing};
    Boundary boundaryCondition;
    positionType treatBoundaries(positionType pos, int individual);
    // inline for speed, may be negative:
    inline positionType getDistance(positionType pos1, positionType pos2) {
        positionType dx = pos1-pos2;
        if (boundaryCondition==Circular) {
            if (dx>maxPos/2) {
                return maxPos-dx;
            } else if (dx<-maxPos/2) {
                return maxPos+dx;
            } else
                return dx;
        } else {
            return dx;
        }
    }
public:
    Continuous_space(Population& p, ParameterFile& pf);
    virtual ~Continuous_space();
    virtual double getPosition(int individual, int dim);
    //positionType& getCoord(int individual, int dim); // position with reference
    virtual void initialize();
    virtual void disperse();
    virtual void prepareNewGeneration(int size);
    virtual void nextGeneration();
    virtual void addChild(int mom, int dad);
    virtual void compactData(std::vector<bool>& alive);
    virtual void resumeFromCheckpoint(Space_sample* S_sample);
    virtual double getDist2(int ind1, int ind2); // squared distance
    virtual double get_mean(int dim);
    virtual double get_variance(int dim, double mean);
    virtual bool isDiscrete();
    virtual Space_sample* make_sample();
};

class Continuous_space_sample : public Space_sample {
    friend class Continuous_space;
    protected:
    Matrix<Continuous_space::positionType> positions;
    public:
    Continuous_space_sample(Continuous_space& S);
    Continuous_space_sample(iSimfile& isf);
    virtual ~Continuous_space_sample();
    virtual void write_to_file(oSimfile& osf);
};
#endif /* defined(__Species__space__) */
