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
#include <assert.h>
#include <iostream>
#include <cmath>
#include <cstring>

#include "mating.h"
#include "fitness.h"
#include "population.h"
#include "randgen.h"
#include "fastmath.h"

// Evolving mate chioces
// 1. f=x, m=x. Magic trait
// 2. f=f, m=x. Evolving Preference for x
// 3. f=y, m=y. Assortment on neutral y
// 4. f=f, m=m. One preference trait, one marker trait



///////////////////////////////////
// Preference
///////////////////////////////////
Preference::Preference(Population& p, ParameterFile & pf) : pop(p)
{
    //Target : Target trait name
    //Preference : Preference trait name (may be same as target)
    //Strength : Strength trait name (may be a constant)

    std::string targetName = pf.getString("display");
    target = pop.findTrait(targetName);
    std::string prefName = pf.getString("preference");
    pref = pop.findTrait(prefName);
    std::string strengthName = pf.getString("strength");
    strength = pop.findTrait(strengthName);
    if ((target->getDims() != pref->getDims()) || (target->getDims() != strength->getDims()) ) {
        std::cout << "Target_selection error: Display, Preference and Strength traits must have the same number of dimensions\n";
        exit(0);
    }
    disassortative_limit = pf.getPositiveDouble("disassortative_limit");
    dislimsq = disassortative_limit*disassortative_limit;
}

double Preference::getPartnerWeight(int female, int male) {
    Fastexp fexp;

    double esum = 0;
    for (int d=0; d<pref->getDims(); ++d) {
        double dx = getY(female,d) - getX(male,d);
        double cd = getC(female,d);
        if (cd>=0) {
            esum += cd * dx * dx;
        } else {
            // disassortative mating
            esum += cd * (dx*dx - dislimsq);
        }
    }
    if (esum>0) {
        return fexp.exp(-esum);
    } else {
        return 1;
    }
}

