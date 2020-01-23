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
#include <cassert>

#include "space.h"
#include "population.h"
#include "randgen.h"


Space::Space(Population& p) : pop(p), Dims(0) {}

Space::~Space() {}

/*class Space_sample {
    protected:
    int pop_size;
    int dims;
    Space_sample(Space& S);
    public:
    virtual void write_to_file(oSimfile& osf)=0;
};*/

Space_sample::Space_sample(Space& S) {
    pop_size = S.pop.size();
    dims = S.getDims();
}

Space_sample::~Space_sample(){
}

void Space_sample::write_to_file(oSimfile &osf) {
    osf.write<size_type>(pop_size);
    osf.write<size_type>(dims);
}

Space_sample::Space_sample(iSimfile& isf) {
    pop_size = isf.read<size_type>();
    dims = isf.read<size_type>();
}

Space_sample* Space_sample::create_from_file(iSimfile &isf) {
    Space_types type = (Space_types)isf.read<char>();
    switch (type) {
        case Space_types::Continuous_space_type:
        return new Continuous_space_sample(isf);
        break;

        case Space_types::Discrete_space_type:
        return new Discrete_space_sample(isf);
        break;

        default:
        std::cout << "Wrong space type!\n";
        exit(1);
        return NULL;
        break;
    }
}

Discrete_space::Discrete_space(Population& p) :
Space(p) {}

bool Discrete_space::isDiscrete() { return true; }

Space_sample* Discrete_space::make_sample() {
    return new Discrete_space_sample(*this);
}

/*class DiscreteSpace_sample {
    protected:
    std::vector<int> linear_patches;
    public:
    DiscreteSpace_sample(DiscreteSpace& S);
}*/

Discrete_space_sample::Discrete_space_sample(Discrete_space& S) : Space_sample(S) {
    linear_patches.clear();
    if (dims>0) { // S may be a NullSpace
        linear_patches.reserve(pop_size);
        for (int ind=0; ind<pop_size; ++ind) {
            linear_patches.push_back(S.getLinearPatch(ind));
        }
    }
}

Discrete_space_sample::~Discrete_space_sample() {
}

void Discrete_space_sample::write_to_file(oSimfile &osf) {
    osf.write<char>((char)Discrete_space_type);
    Space_sample::write_to_file(osf); // writes basic parameters
    osf.writeVector<int>(linear_patches);
}

Discrete_space_sample::Discrete_space_sample(iSimfile& isf) :
Space_sample(isf){
    isf.readVector<int>(linear_patches);
}

Null_space::Null_space(Population& p) : Discrete_space(p) {}
Null_space::~Null_space() {}
void Null_space::initialize() {}
void Null_space::disperse() {}
void Null_space::prepareNewGeneration(int size) {}
void Null_space::nextGeneration() {}
void Null_space::addChild(int mom, int dad) {}
void Null_space::compactData(std::vector<bool>& alive) {}
void Null_space::resumeFromCheckpoint(Space_sample* S_s) {}
double Null_space::getPosition(int individual, int dim) { return 0.0; }
double Null_space::getDist2(int ind1, int ind2) {return 0.0;}

int Null_space::getLinearPatch(int individual) { return 0; }
int Null_space::popSizeInPatch(int linearPatch) {
    if (linearPatch==0) {
        return pop.size();
    } else {
        return 0;
    }
}
int Null_space::getPatchCount() { return 1; }
double Null_space::get_mean(int dim) { return 0.0; }
double Null_space::get_variance(int dim, double mean) { return 0.0; };


Discrete_space_imp::Discrete_space_imp(Population& p, ParameterFile& pf) :
Discrete_space(p) {
    length = pf.getPositiveInt("size");
    Dims = pf.getPositiveInt("dimensions");
    std::string PdTraitname = pf.getString("p_disperse");
    PDisp = pop.findTrait(PdTraitname);
    if (PDisp->get_dims() > 1) {
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

    patches.reserve(Dims,10000);
    newPatches.reserve(Dims,10000);
    patchPopSizes.assign(getPatchCount(), 0);
    generationLastCount = -1;
}

Discrete_space_imp::~Discrete_space_imp() {
}


int Discrete_space_imp::getLinearPatch(int individual) {
  return linearPatches[individual];
}

double Discrete_space_imp::getPosition(int indiviudal, int dim) {
    return patches(dim,indiviudal);
}

double Discrete_space_imp::getDist2(int ind1, int ind2) { // squared distance
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

int Discrete_space_imp::popSizeInPatch(int linearPatch) {
    if (pop.getAge()!=generationLastCount) {
        patchPopSizes.assign(getPatchCount(), 0);
        for (int i=0; i<pop.size(); ++i) {
            patchPopSizes.at(getLinearPatch(i)) += 1;
        }
        generationLastCount = pop.getAge();
    }
    return patchPopSizes.at(linearPatch);
}


void Discrete_space_imp::initialize() {
    // Put everyone in initialPatch
    patches.assign(Dims, pop.size(), initialPatch);
    assignLinearPatches();
}

void Discrete_space_imp::disperse() {
    if (length>1) {
        for (int i=0; i<pop.size(); ++i) {
            if (rand1()<PDisp->traitValue(i)) {
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
                linearPatches[i] = calcLinearPatch(i);
            }
        } // for (int i=0; i<pop.size(); ++i) {
        //assignLinearPatches();
    } // if (length>1)
}

int Discrete_space_imp::treatBoundaries(int pos, int individual) {
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

void Discrete_space_imp::prepareNewGeneration(int size) {
    newPatches.assign(getDims(), 0, 0);
}

void Discrete_space_imp::nextGeneration() {
    patches = newPatches;
    assignLinearPatches();
}

int Discrete_space_imp::calcLinearPatch(int individual) {
    int p = 0;
    for (int d=0; d<Dims; ++d) {
        p = length*p + getPatch(individual, d);
    }
    return p;
}

std::vector<int> Discrete_space_imp::calc_coords_from_linear(int linear) {
    std::vector<int> coords(Dims,0);
    for (int d=Dims-1; d>=0; --d) {
        coords.at(d) = linear % length;
        linear = linear/length;
    }
    return coords;
}

void Discrete_space_imp::assignLinearPatches() {
    if (Dims==1) {
        linearPatches.resize(patches.get_N());
        std::memcpy(&linearPatches.at(0), patches.getX(), patches.get_N()*sizeof(int));
    } else {
        linearPatches.assign(patches.get_N(), 0);
        for (int i=0; i<patches.get_N(); ++i) {
            linearPatches[i] = calcLinearPatch(i);
        }
    }
}

void Discrete_space_imp::addChild(int mom, int dad) {
    // inherit mother's patch:
    newPatches.add_column(&getPatch(mom,0));
}

void Discrete_space_imp::compactData(std::vector<bool>& alive) {
    patches.compact_data(alive);
    int iw=0;
    for (int ir=0; ir<pop.size(); ++ir) {
        if (alive.at(ir)) {
            linearPatches[iw] = linearPatches[ir];
            ++iw;
        }
    }
    linearPatches.resize(iw);
}

void Discrete_space_imp::resumeFromCheckpoint(Space_sample* S_Sample) {
    Discrete_space_sample* DS_sample = dynamic_cast<Discrete_space_sample*>(S_Sample);
    linearPatches = DS_sample->linear_patches;
    assert(linearPatches.size()==pop.size());
    patches.assign(Dims, pop.size(), 0);
    for (int ind=0; ind<pop.size(); ++ind) {
        std::vector<int> coords = calc_coords_from_linear(linearPatches.at(ind));
        for (int d=0; d<Dims; ++d) {
            patches(d,ind) = coords.at(d);
        }
    }
    patchPopSizes.assign(getPatchCount(), 0);
    generationLastCount = -1;
}

double Discrete_space_imp::get_mean(int dim) {
    return patches.row_mean(dim);
}
double Discrete_space_imp::get_variance(int dim, double mean) {
    return patches.row_variance(dim, mean);
}

Continuous_space::Continuous_space(Population& p, ParameterFile& pf) :
Space(p) {
    maxPos = pf.getPositiveDouble("size");
    Dims = pf.getPositiveInt("dimensions");
    std::string PdTraitname = pf.getString("p_disperse");
    PDisp = pop.findTrait(PdTraitname);
    if (PDisp->get_dims() > 1) {
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
    pos.reserve(Dims,10000);
    newPos.reserve(Dims,10000);
}

Continuous_space::~Continuous_space() {
}

bool Continuous_space::isDiscrete() { return false; }

double Continuous_space::getPosition(int indiviudal, int dim) {
    return pos(dim,indiviudal);
}

double Continuous_space::getDist2(int ind1, int ind2) { // squared distance
    if (Dims==1) {
        double dx = pos(0,ind1)-pos(0,ind2);
        return dx*dx;
    } else {
        double dist2 = 0.0;
        positionType *p1 = &pos(0,ind1);
        positionType *p2 = &pos(0,ind2);
        for (int d=0; d<Dims; ++d) {
            double dx = getDistance(*p1, *p2);
            dist2 += dx*dx;
            ++p1;
            ++p2;
        }
        return dist2;
    }
}

void Continuous_space::initialize() {
    // Put everyone in initialPosition
    pos.assign(Dims,pop.size(), initialPosition);
}

Continuous_space::positionType Continuous_space::treatBoundaries(Continuous_space::positionType pos, int individual) {
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

void Continuous_space::disperse() {
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
                pos(0,i) = treatBoundaries(pos(0,i) + dist, i);
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
                    pos(d,i) = treatBoundaries(pos(d,i) + dir[d]*C, i);
                }
            }
        }
    }
}

void Continuous_space::prepareNewGeneration(int size) {
    newPos.assign(getDims(), 0, 0);
}

void Continuous_space::nextGeneration() {
    pos = newPos;
}

void Continuous_space::addChild(int mom, int dad) {
    // inherit mother's position:
    newPos.add_column(&pos(0,mom));
}

void Continuous_space::compactData(std::vector<bool>& alive) {
    pos.compact_data(alive);
}

void Continuous_space::resumeFromCheckpoint(Space_sample* S_sample) {
    Continuous_space_sample* CS_sample = dynamic_cast<Continuous_space_sample*>(S_sample);
    pos = CS_sample->positions;
}

double Continuous_space::get_mean(int dim) {
    return pos.row_mean(dim);
}
double Continuous_space::get_variance(int dim, double mean) {
    return pos.row_variance(dim, mean);
}

Continuous_space_sample::Continuous_space_sample(Continuous_space& S) : Space_sample(S) {
    positions = S.pos;
}

void Continuous_space_sample::write_to_file(oSimfile& osf) {
    osf.write<char>((char)Continuous_space_type);
    Space_sample::write_to_file(osf);
    assert(positions.get_M()==dims);
    assert(positions.get_N()==pop_size);
    positions.write_to_file(osf);
}

Continuous_space_sample::Continuous_space_sample(iSimfile& isf) :
Space_sample(isf){
    positions.read_from_file(isf);
}

Continuous_space_sample::~Continuous_space_sample() {
}

Space_sample* Continuous_space::make_sample() {
    return new Continuous_space_sample(*this);
}
