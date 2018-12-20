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


#include <math.h>
#include <iostream>
#include <algorithm>

#include "fitness.h"
#include "population.h"
#include "fastmath.h"
#include "randgen.h"

Fitness::Fitness(Population& p) : pop(p) {
    //fitness.clear();
    //fitness.reserve(10000);
}

Fitness::~Fitness() {}


StabilizingSelection::StabilizingSelection(Population& p, ParameterFile& pf) :
Fitness(p) {
    std::string zTraitname = pf.getString("trait");
    zTrait = pop.findTrait(zTraitname);
    cost_coefficient = pf.getDouble("cost_coefficient");
    cost_exponent = pf.getDouble("cost_exponent");
}

StabilizingSelection::~StabilizingSelection() {}

void StabilizingSelection::aggregateFitness(std::vector<double>& fitness) {
    //double* fi = &fitness[0];
    for (int i=0; i<pop.size(); ++i) {
        if (fitness[i]>0) {
            double dzsum = 0;
            for (int d=0; d<zTrait->getDims(); ++d) {
                double z = zTrait->traitValue(i,d);
                dzsum += z*z;
            }
            if(cost_exponent!=2) {
                dzsum = fastabspow(dzsum,cost_exponent/2);
            }
            fitness[i] *= std::max(0.0,(1-cost_coefficient*dzsum));
        }
    }
}

DensityDependence::DensityDependence(Population& pop, ParameterFile& pf) :
Fitness(pop) {
    r = pf.getPositiveDouble("r");
    K = pf.getPositiveDouble("k");
    s_space = pf.getDouble("s_space");
    if (s_space<0 ) {
        std::cout << "Error! \n" << "Fitness module Density_dependence: \n" << "Parameter s_space can not be negative.\n";
        exit(0);
    }
    if (s_space==0 && !pop.getSpace().isDiscrete()) {
        std::cout << "Error! \n" << "Fitness module Density_dependence: \n" << "Parameter s_space can not be zero in continuous space.\n";
        exit(0);
    }
}

DensityDependence::~DensityDependence() {}

void DensityDependence::aggregateFitness(std::vector<double>& fitness) {

    Fastexp fexp; // Object for fast exponential function

    // Adjust density dependece for F:
    double KF = K*pop.getF();
    
    // NullSpace case:
    if (pop.getSpace().getDims()==0) {
        int n = pop.size();
        for (int i=0; i<n; ++i) {
            fitness[i] *= fexp.exp( r*(1.0 - n/KF) );
        }
        
    // Within patch case:
    } else if (pop.getSpace().isDiscrete() && s_space==0) {
        DiscreteSpace& theSpace = dynamic_cast<DiscreteSpace&>(pop.getSpace());
        int n = pop.size();
        
        // Vector of local patch densities (linear indexing):
        std::vector<int> N(theSpace.getPatchCount(), 0);
        for (int i=0; i<n; ++i) {
            N[theSpace.getLinearPatch(i)] += 1;
        }
        for (int i=0; i<n; ++i) {
            fitness[i] *= fexp.exp( r*(1.0 - N[theSpace.getLinearPatch(i)]/KF) );
        }
        
    // General case:
    } else {
        double ss2 = s_space*s_space;
        for (int i=0; i<pop.size(); ++i) {
            double Neff = 0; // efficient number of competitors
            for (int ci=0; ci<pop.size(); ++ci) {
                double d2exp = pop.getSpace().getDist2(i, ci) / ss2 / 2;
                if (d2exp < 14) { // exclude contributions smaller than 1e-6
                    Neff += fexp.exp(-d2exp);
                }
            }
            fitness[i] *= fexp.exp( r*(1.0 - Neff/KF) );
        }
    }
}

ResourceLandscape::ResourceLandscape(Population& pop, ParameterFile& pf) :
Fitness(pop)
{
    std::string xTraitname = pf.getString("trait");
    xTrait = pop.findTrait(xTraitname);
    r = pf.getPositiveDouble("r");
    K0 = pf.getPositiveDouble("k_0");
    sK = pf.getPositiveDouble("s_k");
    sa = pf.getPositiveDouble("s_a");
    s_space = pf.getDouble("s_space");
    if (s_space<0 ) {
        std::cout << "Error! \n" << "Fitness module Resource_landscape: \n" << "Parameter s_space can not be negative.\n";
        exit(0);
    }
    if (s_space==0 && !pop.getSpace().isDiscrete()) {
        std::cout << "Error! \n" << "Fitness module Resource_landscape: \n" << "Parameter s_space can not be zero in continuous space.\n";
        exit(0);
    }
    k_space = pf.getDouble("k_space");
}

ResourceLandscape::~ResourceLandscape() {
    //delete [] comp;
}

double ResourceLandscape::getTraitDist2( int ind1, int ind2) { // squared distance in trait space
    double dx2 = 0;
    int dims = xTrait->getDims();
    if (dims==1) {
        dx2 = getX(ind1) - getX(ind2);
        dx2 *= dx2;
    } else {
        double *xp1 = &getX(ind1,0);
        double *xp2 = &getX(ind2,0);
        for (int d=0; d<dims; ++d) {
            double dx = *xp1 - *xp2;
            dx2 += dx*dx;
            ++xp1;
            ++xp2;
        }
    }
    return dx2;
}


void ResourceLandscape::aggregateFitness(std::vector<double>& fitness)
 {
     double twosa2 = 2*sa*sa;
     int n = pop.size();
     int xdims = xTrait->getDims();
     int sdims = pop.getSpace().getDims();

     std::vector<double> comp;
     comp.assign(n, 0); // vector of efficient competition, one per individual
     Fastexp fexp;
     
     // Calculate effective competition
     if (xdims==1 && sdims==0) {
         // Special case 1: 1D trait, no space
         for (int i=0; i<n; ++i) {
             double compi = comp.at(i);
             // self competition:
             compi += 1;
             if (i<n-1) {
                 double ix = getX(i);
                 // This is a loop of complexity n^2. Use efficient iterators.
                 double* xcip = &getX(i+1);
                 double* compcip = &comp.at(i+1);
                 for (int ci=i+1; ci<n; ++ci) {
                     double xdiff = ix - *xcip;
                     double c = fexp.exp(-xdiff*xdiff/twosa2 );
                     compi += c;
                     *compcip += c;
                     ++compcip;
                     ++xcip;
                 }
             }
             comp.at(i) = compi;
         }
         for (int i=0; i<n; ++i) {
             double ix = getX(i);
             double iKF = K0*fexp.exp(-ix*ix/2/sK/sK) * pop.getF(); // adjustment such that K is carrying capacity AFTER selection
             fitness.at(i) *= std::max(0.0,1 + r*(1-comp.at(i)/iKF));
         }
         return;
     } else if (pop.getSpace().isDiscrete() && s_space==0) {
         // Special case 2: within patch competition
         DiscreteSpace& theSpace = dynamic_cast<DiscreteSpace&>(pop.getSpace());
         for (int i=0; i<n; ++i) {
             double compi = comp.at(i);
             int patchi = theSpace.getLinearPatch(i);
             // self competition:
             compi += 1;
             if (i<n-1) {
                 double* compcip = &comp.at(i+1);// pointer to competitor's comp-value
                 for (int ci=i+1; ci<n; ++ci) {
                     if (theSpace.getLinearPatch(ci) == patchi) {
                         double x2sum = getTraitDist2(i, ci);
                         double c = fexp.exp(-x2sum/twosa2 );
                         compi += c;
                         *compcip += c;
                     }
                     ++compcip;
                 }
             }
             comp.at(i) = compi;
         }
     }
     else { // general case
         double twoss2 = 2*s_space*s_space;
         for (int i=0; i<n; ++i) {
             double compi = comp.at(i);
             // self competition:
             compi += 1;
             if (i<n-1) {
                 for (int ci=i+1; ci<n; ++ci) {
                     double dist2 = pop.getSpace().getDist2(i, ci);
                     double x2sum = getTraitDist2(i, ci);
                     double c = fexp.exp(-x2sum/twosa2 - dist2/twoss2);
                     compi += c;
                     comp.at(ci) += c;
                 }
             }
             comp.at(i) = compi;
         }
     }
     if (xdims==1) {
         double *xip = &getX(0,0);
         for (int i=0; i<n; ++i) {
             double xopt = k_space*pop.getSpace().getPosition(i, 0);
             double dx = *xip - xopt;
             double iKF = K0*fexp.exp(-dx*dx/2/sK/sK) * pop.getF(); // adjustment such that K is carrying capacity AFTER selection
             fitness.at(i) *= std::max(0.0,1 + r*(1-comp.at(i)/iKF));
             ++xip;
         }
     } else {
         double *xip = &getX(0,0);
         for (int i=0; i<n; ++i) {
             double xopt = k_space*pop.getSpace().getPosition(i, 0);
             double dx2 = 0;
             for (int d=0; d<xdims; ++d) {
                 double dx = *xip - xopt;
                 dx2 += dx*dx;
                 ++xip;
             }
             //double ix = getX(i) - k_space*pop.getSpace().getPosition(i, 0);
             double iKF = K0*fexp.exp(-dx2/2/sK/sK) * pop.getF(); // adjustment such that K is carrying capacity AFTER selection
             fitness.at(i) *= std::max(0.0,1 + r*(1-comp.at(i)/iKF));
         }
     }
}


DiscreteResources::DiscreteResources(Population& pop, ParameterFile& pf) :
Fitness(pop)
{
    if (!pop.getSpace().isDiscrete()) {
        std::cout << "Dicrete_resources module is not implemented for continuous space.\n";
        exit(0);
    }
    std::string xTraitname = pf.getString("trait");
    xTrait = pop.findTrait(xTraitname);
    nR = pf.getPositiveInt("n_r");
    K = pf.getPositiveDouble("k");
    a0 = pf.getPositiveDouble("a_0");
    double sa = pf.getPositiveDouble("s_a");
    ta = 1/sa/sa;
    cmin = pf.getDouble("c_min");
}

void DiscreteResources::aggregateFitness(std::vector<double>& fitness) {
    /* model:
     Ri = K/(1+sum(aij*Nj/K))
     fitness_j = 1 + sum_i(aij*Ri) - cmin;
     */
    Fastexp fexp;
    
    // attack rates:
    std::vector<double> a(pop.size());
    
    // Individual intakes:
    std::vector<double> intake(pop.size());
    intake.assign(pop.size(), 0);
    
    if (pop.getSpace().getDims()==0) { // no space
        for (int ri=0; ri<nR; ++ri) {
            double consumption = 0;
            for (int i=0; i<pop.size(); ++i) {
                double dx = getX(i) - ri;
                a[i] = a0 * fexp.exp(-ta*dx*dx/2.0) / K;
                consumption += a[i]; // should be divided by F. Adjusted outside loop.
            }
            consumption /= pop.getF(); // adjust consumer density for F.
            // Resource abundance:
            double R = K/(1 + consumption);
            for (int i=0; i<pop.size(); ++i) {
                intake[i] += a[i]*R;
            }
        }
    } else { // we can assume discrete space
        DiscreteSpace& theSpace = dynamic_cast<DiscreteSpace&>(pop.getSpace());
        for (int patchi=0; patchi<theSpace.getPatchCount(); ++patchi) {
            for (int ri=0; ri<nR; ++ri) {
                double consumption = 0;
                for (int i=0; i<pop.size(); ++i) {
                    if( theSpace.getLinearPatch(i) == patchi ) {
                        double dx = getX(i) - ri;
                        a[i] = a0*fexp.exp(-ta*dx*dx/2.0) / K;
                        consumption += a[i]; // should be divided by F. Adjusted outside loop.
                    }
                }
                consumption /= pop.getF(); // adjust consumer density for F.
                // Resource abundance:
                double R = K/(1 + consumption);
                for (int i=0; i<pop.size(); ++i) {
                    if( theSpace.getLinearPatch(i) == patchi ) {
                        intake[i] += a[i]*R;
                    }
                }
            }
        }
    }

    for (int i=0; i<pop.size(); ++i) {
        fitness[i] *= std::max(0.0, 1 + intake[i] - cmin);
    }
}

DiscreteResources::~DiscreteResources() {
}


////////////////////////////////// Spatial Gradient

SpatialGradient::SpatialGradient(Population& pop, ParameterFile& pf) :
Fitness(pop) {
    std::string xTraitname = pf.getString("trait");
    xTrait = pop.findTrait(xTraitname);
    ks = pf.getDouble("k_space");
    double ss = pf.getPositiveDouble("s_selection");
    ts = 1/ss/ss;
}

SpatialGradient::~SpatialGradient() {}

void SpatialGradient::aggregateFitness(std::vector<double>& fitness) {
    Fastexp fexp;
    for (int i=0; i<pop.size(); ++i) {
        // adaptation in first dimension only:
        double dx = xTrait->traitValue(i,0) - ks*pop.getSpace().getPosition(i,0);
        fitness[i] *= fexp.exp(-ts*dx*dx/2.0);
    }
}

Catastrophe::Catastrophe(Population& pop, ParameterFile& pf) :
Fitness(pop) {
    Pcat = pf.getDouble("p_catastrophe");
    if (Pcat<=0) {
        std::cout << "Warning, Catastrophe module: there will be no catastrophes. P_catastrophe = " << Pcat << '\n';
    }
    Psurv = pf.getPositiveDouble("p_survive");
    if (Psurv<=0) {
        std::cout << "Warning, Catastrophe module: Catastrophes will cause global extinction. P_survive = " << Psurv << '\n';
    }
}

Catastrophe::~Catastrophe() {}

void Catastrophe::aggregateFitness(std::vector<double>& fitness) {
    if (rand1()<Pcat) {
        for (int i=0; i<pop.size(); ++i) {
            fitness[i] *= Psurv;
        }
    }
}

