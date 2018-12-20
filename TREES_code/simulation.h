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


#ifndef __Species__simulation__
#define __Species__simulation__

#include <stdio.h>
#include <string>
#include <vector>
#include <sstream>

#include "population.h"

/////////////////////////////
// Simulation class
// The Simulation class is the overarching keeper of simulation parameters
// and simulation output (samples).
// It runs the simulation through a simple for-loop over time and
// eventually saves all parameters and population samples to a 'sim'-file.
///////////////////////////////////////////////
class Simulation {
protected:
    unsigned seed;
    bool verbose;
    timeType tmax;
    timeType sampleInterval;
    std::string  simName;
    std::ostringstream  parameterFileCopy;
    Population* thePop;
    std::vector<Sample*> sampleList;
    void save(int replicate);

public:
    Simulation(std::string& parameterFileName);
    ~Simulation();
    void runAndSave(int iterations);
};
#endif /* defined(__Species__simulation__) */
