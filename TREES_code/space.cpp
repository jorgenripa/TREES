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


#include <cctype>
#include <iostream>
#include <math.h>
#include "stdlib.h"

#include "space.h"
#include "population.h"
#include "randgen.h"


Space::Space(Population& p) : pop(p), Dims(0) {}

Space::~Space() {}

DiscreteSpace::DiscreteSpace(Population& p) :
Space(p) {}

bool DiscreteSpace::isDiscrete() { return true; }

NullSpace::NullSpace(Population& p) : DiscreteSpace(p) {}
NullSpace::~NullSpace() {}
void NullSpace::initialize(int n0) {}
void NullSpace::disperse() {}
void NullSpace::prepareNewGeneration(int size) {}
void NullSpace::nextGeneration() {}
void NullSpace::addChild(int mom, int dad) {}
void NullSpace::compactData(std::vector<bool>& alive) {}
void NullSpace::addToSample(Sample& theSample) {}
double NullSpace::getPosition(int individual, int dim) { return 0.0; }
double NullSpace::getDist2(int ind1, int ind2) {return 0.0;}

int NullSpace::getLinearPatch(int individual) { return 0; }
int NullSpace::popSizeInPatch(int linearPatch) {
    if (linearPatch==0) {
        return pop.size();
    } else {
        return 0;
    }
}
int NullSpace::getPatchCount() { return 1; }


DiscreteSpaceImp::DiscreteSpaceImp(Population& p, ParameterFile& pf) :
DiscreteSpace(p) {
    length = pf.getPositiveInt("size");
    Dims = pf.getPositiveInt("dimensions");
    std::string PdTraitname = pf.getString("p_disperse");
    PDisp = pop.findTrait(PdTraitname);
    if (PDisp->getDims() > 1) {
        std::cout << "Dispersal probability can not be multidimensional. Trait : " << PdTraitname << '\n';
        exit(0);
    }
    std::string dispS = pf.getString("dispersal_type");
    switch (std::toupper(dispS[0])) {
        case 'G': dispType = Global; break;
        case 'N': dispType = Neighbour; break;
        case 'D':
            dispType = Distance;
            dispdist = pf.getPositiveDouble("dispersal_distance");
            Pjump = 1-1/dispdist;
            break;
        default:
            std::cout << "Unknown dispersal_type : " << dispS << '\n';
            exit(0);
            break;
    }
    
    std::string bounds = pf.getString("boundary");
    switch (std::toupper(bounds[0])) {
        case 'R': boundaryCondition = Reflective; break;
        case 'C': boundaryCondition = Circular; break;
        case 'A': boundaryCondition = Absorbing; break;
        default:
            std::cout << "Unknown boundary condition : " << bounds << '\n';
            exit(0);
    }
    
    initialPatch = pf.getInt("initial_position");
    if (initialPatch<0 || initialPatch>=length) {
        std::cout << "DiscreteSpace: Initial position outside range : " << initialPatch << '\n';
        exit(0);
    }

    patches.reserve(10000*Dims);
    newPatches.reserve(10000*Dims);
    patchPopSizes.assign(getPatchCount(), 0);
    generationLastCount = -1;
}

DiscreteSpaceImp::~DiscreteSpaceImp() {
}


int DiscreteSpaceImp::getLinearPatch(int individual) {
  return linearPatches[individual];
}

double DiscreteSpaceImp::getPosition(int indiviudal, int dim) {
    return patches[indiviudal*Dims + dim];
}

double DiscreteSpaceImp::getDist2(int ind1, int ind2) { // squared distance
    double dist2=0.0;
    int *p1 = &getPatch(ind1, 0);
    int *p2 = &getPatch(ind2, 0);
    for (int d=0; d<Dims; ++d) {
        int dx = getDistance(*p1, *p2);
        dist2 += dx*dx;
        ++p1;
        ++p2;
    }
    return dist2;
}

int DiscreteSpaceImp::popSizeInPatch(int linearPatch) {
    if ((int)pop.getAge()!=generationLastCount) {
        patchPopSizes.assign(getPatchCount(), 0);
        for (int i=0; i<pop.size(); ++i) {
            patchPopSizes.at(getLinearPatch(i)) += 1;
        }
        generationLastCount = (int)pop.getAge();
    }
    return patchPopSizes.at(linearPatch);
}


void DiscreteSpaceImp::initialize(int n0) {
    // Put everyone in initialPatch
    patches.assign(n0*Dims, initialPatch);
    if (Dims==1) {
        linearPatches = patches;
    } else {
        int p = 0;
        for (int d=0; d<Dims; ++d) {
            p = length*p + initialPatch;
        }
        linearPatches.assign(patches.size(), p);
    }
}

void DiscreteSpaceImp::disperse() {
    if (length>1) {
        for (int i=0; i<pop.size(); ++i) {
            if (rand1()<PDisp->traitValue(i) ) {
                double z;
                int dim, dist;
                switch (dispType) {
                    case Global:
                        for (int d=0; d<Dims; ++d) {
                            getPatch(i, d) = rand1()*(length-1);
                        }
                        break;
                    case Neighbour:
                        z = rand1()*Dims;
                        // choose one random dimension for dispersal:
                        dim = int(z);
                        // recycle z to choose direction:
                        if (z-dim<0.5) {
                            getPatch(i, dim) = treatBoundaries(getPatch(i, dim)+1, i);
                        } else {
                            getPatch(i, dim) = treatBoundaries(getPatch(i, dim)-1, i);
                        }
                        break;
                    case Distance:
                        dist=1;
                        // Dispersal distance has a geometric distribution.
                        // Mean distance is dispdist.
                        while (rand1()<Pjump) {
                            dist += 1;
                        }
                        z = rand1()*Dims;
                        // choose one random dimension for dispersal:
                        dim = int(z);
                        // recycle z to choose direction:
                        if (z-dim<0.5) {
                            getPatch(i, dim) = treatBoundaries(getPatch(i, dim)+dist, i);
                        } else {
                            getPatch(i, dim) = treatBoundaries(getPatch(i, dim)-dist, i);
                        }
                        break;
                    default:
                        break;
                }
            }
        }
    }
}

int DiscreteSpaceImp::treatBoundaries(int pos, int individual) {
    if (pos>=length || pos<0) {
        switch (boundaryCondition) {
            case Reflective:
                // reflect until within bounds:
                while (pos<0 || pos>=length) {
                    if (pos<0) {
                        pos = -pos;
                    } else {
                        pos = 2*(length-1)-pos;
                    }
                }
                break;
            case Circular:
                pos = pos % length;
                if (pos<0) {
                    pos += length;
                }
                break;
            case Absorbing:
                pop.kill(individual);
                break;
            default:
                // this should never happen
                break;
        }
    }
    return pos;
}


void DiscreteSpaceImp::prepareNewGeneration(int size) {
    newPatches.clear();
}

void DiscreteSpaceImp::nextGeneration() {
    patches = newPatches;
    if (Dims==1) {
        linearPatches = patches;
    } else {
        linearPatches.assign(patches.size(), 0);
        for (int i=0; i<patches.size(); ++i) {
            int p = 0;
            for (int d=0; d<Dims; ++d) {
                p = length*p + getPatch(i, d);
            }
            linearPatches[i] = p;
        }
    }
}

void DiscreteSpaceImp::addChild(int mom, int dad) {
    // inherit mother's patch:
    for (int d=0; d<Dims; ++d) {
        newPatches.push_back(getPatch(mom,d));
    }
}

void DiscreteSpaceImp::compactData(std::vector<bool>& alive) {
    int iw=0;
    for (int ir=0; ir<pop.size(); ++ir) {
        if (alive.at(ir)) {
            for (int d=0; d<Dims; ++d) {
                getPatch(iw, d) = getPatch(ir, d);
            }
            linearPatches[iw] = linearPatches[ir];
            ++iw;
        }
    }
    patches.resize(iw*Dims);
    linearPatches.resize(iw);
}

void DiscreteSpaceImp::addToSample(Sample& theSample) {
    // add sample spatial position:
    theSample.addData(new IntData(patches));
}


ContinuousSpace::ContinuousSpace(Population& p, ParameterFile& pf) :
Space(p) {
    maxPos = pf.getPositiveDouble("size");
    Dims = pf.getPositiveInt("dimensions");
    std::string PdTraitname = pf.getString("p_disperse");
    PDisp = pop.findTrait(PdTraitname);
    if (PDisp->getDims() > 1) {
        std::cout << "Dispersal probability can not be multidimensional. Trait : " << PdTraitname << '\n';
        exit(0);
    }
    
    dispdist = pf.getPositiveDouble("dispersal_distance");
    std::string bounds = pf.getString("boundary");
    switch (std::toupper(bounds[0])) {
        case 'R': boundaryCondition = Reflective; break;
        case 'C': boundaryCondition = Circular; break;
        case 'A': boundaryCondition = Absorbing; break;
        default:
            std::cout << "Unknown boundary condition\n";
            exit(0);
    }
    initialPosition = pf.getDouble("initial_position");
    if (initialPosition<0 || initialPosition>maxPos) {
        std::cout << "ContinuousSpace: Initial position outside range : " << initialPosition << '\n';
        exit(0);
    }
    pos.reserve(10000*Dims);
    newPos.reserve(10000*Dims);
}

ContinuousSpace::~ContinuousSpace() {
}

bool ContinuousSpace::isDiscrete() { return false; }

double& ContinuousSpace::getCoord(int indiviudal, int dim) {
    return pos[indiviudal*Dims + dim];
}

double ContinuousSpace::getPosition(int indiviudal, int dim) {
    return pos[indiviudal*Dims + dim];
}

double ContinuousSpace::getDist2(int ind1, int ind2) { // squared distance
    if (Dims==1) {
        double dx = getCoord(ind1, 0)-getCoord(ind2, 0);
        return dx*dx;
    } else {
        double dist2=0.0;
        double *p1 = &getCoord(ind1, 0);
        double *p2 = &getCoord(ind2, 0);
        for (int d=0; d<Dims; ++d) {
            double dx = getDistance(*p1, *p2);
            dist2 += dx*dx;
            ++p1;
            ++p2;
        }
        return dist2;
    }
}

// this is inline for speed:
/*double ContinuousSpace::getDistance(double pos1, double pos2) {
    double dx = fabs(pos1-pos2);
    if (boundaryCondition==circular && dx>maxPos/2) {
        return maxPos-dx;
    } else {
        return dx;
    }
}*/

void ContinuousSpace::initialize(int n0) {
    // Put everyone in initialPosition
    pos.assign(n0*Dims, initialPosition);
}

double ContinuousSpace::treatBoundaries(double pos, int individual) {
    if (pos>maxPos || pos<0) {
        switch (boundaryCondition) {
            case Reflective:
                while (pos<0 || pos>maxPos) {
                    if (pos<0) {
                        pos = -pos;
                    } else {
                        pos = 2*maxPos-pos;
                    }
                }
                break;
            case Circular:
                if (pos>maxPos) {
                    pos = fmod(pos,maxPos);
                } else {
                    pos = fmod(pos,maxPos)+maxPos;
                }
                break;
            case Absorbing:
                pop.kill(individual);
                break;
            default:
                // this should never happen
                break;
        }
    }
    return pos;
}

void ContinuousSpace::disperse() {
    for (int i=0; i<pop.size(); ++i) {
        double z = rand1();
        double PD = PDisp->traitValue(i);
        if (z<PD) {
            // generate an expontentially distributed dispersal distance:
            double dist = -dispdist*log(rand1());
            if (Dims==1) {
                // recycling of random number:
                if (z<PD/2) {
                    dist = -dist;
                }
                getCoord(i,0) = treatBoundaries(getCoord(i,0) + dist, i);
            } else {
                // Generate a random direction:
                std::vector<double> dir(Dims,0);
                double r2 = 0;
                for (int d=0; d<Dims; ++d) {
                    dir[d] = randn();
                    r2 += dir[d]*dir[d];
                }
                double C = dist/sqrt(r2);
                for (int d=0; d<Dims; ++d) {
                    getCoord(i,d) = treatBoundaries(getCoord(i,d) + dir[d]*C, i);
                }
            }
        }
    }
}

void ContinuousSpace::prepareNewGeneration(int size) {
    newPos.clear();
}

void ContinuousSpace::nextGeneration() {
    pos = newPos;
}

void ContinuousSpace::addChild(int mom, int dad) {
    // inherit mother's position:
    for (int d=0; d<Dims; ++d) {
        newPos.push_back(getCoord(mom,d));
    }
}

void ContinuousSpace::compactData(std::vector<bool>& alive) {
    int iw=0; // writing position
    for (int ir=0; ir<pop.size(); ++ir) {
        if (alive.at(ir)) {
            for (int d=0; d<Dims; ++d) {
                getCoord(iw,d) = getCoord(ir,d);
            }
            ++iw;
        }
    }
    pos.resize(iw*Dims);
}

void ContinuousSpace::addToSample(Sample& theSample) {
    // add sample spatial position:
    theSample.addData(new FloatData(pos));
}

