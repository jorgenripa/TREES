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


#ifndef __Species__fitness__
#define __Species__fitness__

#include <stdio.h>
#include <vector>

#include "trait.h"
#include "parameterFile.h"

/////////////////////////////
// Fitness modules
// Used in combination, multiplicative fitness components
/////////////////////////////

class Population;

// The Fitness interface:
class Fitness {
protected:
    Population& pop;
    Fitness(Population& p);
public:
    virtual ~Fitness();
    virtual void aggregateFitness(std::vector<double>& fitness)=0;
    //virtual double getFitness(int i);
};

class StabilizingSelection : public Fitness {
protected:
    // cost = costCoefficient*sum(|dx|^costExponent)
    Trait* zTrait;
    double optimum;
    double cost_coefficient;
    double cost_exponent;
public:
    StabilizingSelection(Population& p, ParameterFile& pf);
    virtual ~StabilizingSelection();
    virtual void aggregateFitness(std::vector<double>& fitness);
    //virtual double getFitness(int i);
};

class DensityDependence : public Fitness {
protected:
    double r;
    double K;
    double s_space;
public:
    DensityDependence(Population& pop, ParameterFile& pf);
    virtual ~DensityDependence();
    virtual void aggregateFitness(std::vector<double>& fitness);
};


class ResourceLandscape : public Fitness {
protected:
    double r, K0, sK, sa, s_space, k_space;
    //    inline double& getTrait(int individual) { return traits[0].traitValue(individual); }
    Trait* xTrait;
    inline traitType& getX(int individual) { return xTrait->traitValue(individual); }
    inline traitType& getX(int individual, int dim) { return xTrait->traitValue(individual,dim); }
    double getTraitDist2( int ind1, int ind2); // squared distance in trait space
public:
    ResourceLandscape(Population& pop, ParameterFile& pf);
    virtual ~ResourceLandscape();
    virtual void aggregateFitness(std::vector<double>& fitness);
};


// discrete resource model:
class DiscreteResources : public Fitness {
    /* model:
     Ri = K/(1+sum(aij*Nj/K))
     fitness_j = 1 + sum_i(aij*Ri)/K - cmin;
     */
protected:
    int nR; // number of resources
    double K; // system scale
    double a0; // maximal attack rate
    double ta; // trade-off
    double cmin; // minimal consumption level for status quo (standard = 1)
    Trait* xTrait; // resource adaptation trait (only first dimension is used)
    inline traitType& getX(int individual) { return xTrait->traitValue(individual,0); }
public:
    DiscreteResources(Population& pop, ParameterFile& pf);
    virtual ~DiscreteResources();
    virtual void aggregateFitness(std::vector<double>& fitness);
};


// Local adaptation on a spatial gradient (discrete or continuous)
// Only first dimension is used
class SpatialGradient : public Fitness {
protected:
    Trait* xTrait; // local adaptation trait (only first dimension is used)
    double ks;
    double ts; // spatial trade-off strength
public:
    SpatialGradient(Population& pop, ParameterFile& pf);
    virtual ~SpatialGradient();
    virtual void aggregateFitness(std::vector<double>& fitness);
};


class Catastrophe : public Fitness {
protected:
    double Pcat; // probability of catastrophic event, per generation
    double Psurv; // severity of catastrophe, survival probability
public:
    Catastrophe(Population& pop, ParameterFile& pf);
    virtual ~Catastrophe();
    virtual void aggregateFitness(std::vector<double>& fitness);
};

#endif /* defined(__Species__fitness__) */
