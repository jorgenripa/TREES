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
#include "sample.h"
#include "trait.h"

/////////////////////////////
// Space modules
// A Space object holds the spatial as well as dispersal parameters.
// I keeps track of the positions of all individuals and manages dispersal.
/////////////////////////////

class Population;

// The space interface:
class Space {
protected:
    Population& pop;
    int Dims;
public:
    Space(Population& p);
    virtual ~Space();
    int getDims() { return Dims; }
    virtual void initialize(int n0)=0;
    virtual void disperse()=0;
    virtual void prepareNewGeneration(int size)=0;
    virtual void nextGeneration()=0;
    virtual void addChild(int mom, int dad)=0;
    virtual void compactData(std::vector<bool>& alive)=0;
    virtual void addToSample(Sample& theSample)=0;
    virtual double getPosition(int individual, int dimension)=0;
    virtual double getDist2(int ind1, int ind2)=0; // squared distance
    virtual bool isDiscrete()=0;
};

// The DiscreteSpace interface
class DiscreteSpace : public Space {
public:
    DiscreteSpace(Population& p);
    //virtual ~DiscreteSpace()=0;
    // new functions:
    virtual int getLinearPatch(int individual)=0;
    virtual int popSizeInPatch(int linearPatch)=0;
    virtual int getPatchCount()=0;
    // From Space interface:
    virtual void initialize(int n0)=0;
    virtual void disperse()=0;
    virtual void prepareNewGeneration(int size)=0;
    virtual void nextGeneration()=0;
    virtual void addChild(int mom, int dad)=0;
    virtual void compactData(std::vector<bool>& alive)=0;
    virtual void addToSample(Sample& theSample)=0;
    virtual double getPosition(int individual, int dim)=0;
    virtual double getDist2(int ind1, int ind2)=0; // squared distance
    virtual bool isDiscrete();
};

// A discrete space implementation, with global dispersal.
class DiscreteSpaceImp : public DiscreteSpace {
protected:
    std::vector<int> patches;
    std::vector<int> newPatches;
    std::vector<int> patchPopSizes;
    std::vector<int> linearPatches;
    timeType generationLastCount;
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
public:
    DiscreteSpaceImp(Population& p, ParameterFile& pf);
    virtual ~DiscreteSpaceImp();
    int& getPatch(int individual, int dim) {return patches.at(individual*Dims + dim);}
    // DiscreteSpace interface:
    virtual int getLinearPatch(int individual);
    virtual int popSizeInPatch(int linearPatch);
    virtual int getPatchCount() { return std::pow(length,Dims); }
    // Space interface:
    virtual void initialize(int n0);
    virtual void disperse();
    virtual void prepareNewGeneration(int size);
    virtual void nextGeneration();
    virtual void addChild(int mom, int dad);
    virtual void compactData(std::vector<bool>& alive);
    virtual void addToSample(Sample& theSample);
    virtual double getPosition(int individual, int dim);
    virtual double getDist2(int ind1, int ind2); // squared distance
};

// The NullSpace class is a DiscreteSpace implementation
//   with no stored positions.
// All individuals are regarded as positioned in patch 0.
class NullSpace : public DiscreteSpace {
public:
    NullSpace(Population& p);
    virtual ~NullSpace();

    virtual int getLinearPatch(int individual);
    virtual int popSizeInPatch(int linearPatch);
    virtual int getPatchCount();
    // From Space interface:
    virtual void initialize(int n0);
    virtual void disperse();
    virtual void prepareNewGeneration(int size);
    virtual void nextGeneration();
    virtual void addChild(int mom, int dad);
    virtual void compactData(std::vector<bool>& alive);
    virtual void addToSample(Sample& theSample);
    virtual double getPosition(int individual, int dimension);
    virtual double getDist2(int ind1, int ind2); // squared distance
};

/*class CSIterator {
public:
    double* ppos;
    int dims;
    CSIterator& operator ++() { ppos += dims; return *this; }
};*/

class ContinuousSpace : public Space {
protected:
    std::vector<double> pos;
    std::vector<double> newPos;
    double initialPosition;  // starting position for all individuals (in all dimensions)
    Trait* PDisp;
    double dispdist; // mean dispersal distance
    double maxPos; // size in each dimension
    enum Boundary {Circular, Reflective, Absorbing};
    Boundary boundaryCondition;
    double treatBoundaries(double pos, int individual);
    // inline for speed, may be negative:
    inline double getDistance(double pos1, double pos2) {
        double dx = pos1-pos2;
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
    ContinuousSpace(Population& p, ParameterFile& pf);
    virtual ~ContinuousSpace();
    virtual double getPosition(int individual, int dim);
    double& getCoord(int individual, int dim); // position with reference
    virtual void initialize(int n0);
    virtual void disperse();
    virtual void prepareNewGeneration(int size);
    virtual void nextGeneration();
    virtual void addChild(int mom, int dad);
    virtual void compactData(std::vector<bool>& alive);
    virtual void addToSample(Sample& theSample);
    virtual double getDist2(int ind1, int ind2); // squared distance
    virtual bool isDiscrete();
};


#endif /* defined(__Species__space__) */
