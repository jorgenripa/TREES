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


#ifndef __Species__mating__
#define __Species__mating__

#include <stdio.h>
#include <vector>

#include "trait.h"
#include "parameterFile.h"

class Population;
class Survival;

enum MatingType {Selfing, Local, Random};

// Evolving mate chioces
// 1. f=x, m=x. Magic trait
// 2. f=f, m=x. Evolving Preference for x
// 3. f=y, m=y. Assortment on neutral y
// 4. f=f, m=m. One preference trait, one marker trait

//Mating : Preference
// Class with three traits:
//Target : Target trait name
//Preference : Preference trait name (may be same as target)
//Strength : Strength trait name (may be a constant)
// 1. a preferred x (value of arbitrary trait)
// 2. the degree of choosiness
class Preference {
protected:
    Population &pop;
    Trait* target;// X
    Trait* pref; // Y
    Trait* strength; // C
    double disassortative_limit;
    double dislimsq;
    double& getX(int individual,int dim) { return target->traitValue(individual,dim);}
    double& getY(int individual,int dim) { return pref->traitValue(individual,dim);}
    double& getC(int individual, int dim) { return strength->traitValue(individual,dim);}
public:
    Preference(Population& p, ParameterFile & pf);
    double getPartnerWeight(int female, int male);
};


#endif /* defined(__Species__mating__) */
