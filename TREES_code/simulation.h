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
// Checkpoints, if activated, are saved to the same file.
///////////////////////////////////////////////
class Simulation {
protected:
    seed_type seed;
    bool verbose;
    time_type tmax;
    char microsample_option; // N[one], M[eans], V[ariance], C[ovariance] (all converted to lower case)
    time_type sampleInterval;
    time_type checkpointInterval; // how often to save checkpoints
    bool keepOldCheckpoints;
    std::string  simName;
    std::string  resultsFileName;
    std::ostringstream  parameterFileCopy;
    Population* thePop;
    std::vector<Sample_base*> sample_list;
    std::vector<Sample_base*> microsample_list;
    // CheckpointList is implemented as a list of Sample_base* to be able to administrate it as such. The pointers are typecast to Population_checkpoint* (a subclass of Sample_base) when necessary.
    std::vector<Sample_base*> checkpointList;
    void load(std::string& resultsFile);
    void save();
    void makeFileName(int replicate);
public:
    Simulation(std::string& parameterFileName, bool verbose);
    ~Simulation();
    seed_type get_seed() { return seed;}
    void initialize(int replicate);
    void runFromGeneration(time_type gen, seed_type gen_seed, double time_so_far);
    void resume(std::string resultsFile, time_type gen, std::string newResultsFile);
};

#endif /* defined(__Species__simulation__) */
