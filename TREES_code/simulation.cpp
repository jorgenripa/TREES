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
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>

#include "simulation.h"
#include "parameterFile.h"
#include "population.h"
#include "randgen.h"
#include "simfile.h"
#include "sample.h"
#include "mytime.h"

Simulation::Simulation(std::string& parameterFileName){
    ParameterFile pf(parameterFileName);

    simName = parameterFileName;
    // skip file extension:
    simName.erase(simName.find_last_of('.'));
    std::cout << "name: " << simName << '\n';
    verbose = pf.getBool("verbose");
    tmax = pf.getPositiveLong("t_max");
    sampleInterval = pf.getPositiveLong("sample_interval");

    std::string seedS = pf.getStringLower("seed");
    if (seedS[0]=='r') {
        seed = randomize();
    } else {
        seed = (unsigned)std::atol(seedS.c_str());
        randseed(seed);
    }
    
    //Parse population parameters (including gene tracking/sampling):
    thePop = new Population(pf);
    
    pf.close();

    sampleList.clear();

    // save copy of parameterFile:
    std::ifstream f(parameterFileName);
    if (f) {
        parameterFileCopy << f.rdbuf();
        f.close();
    }
}

Simulation::~Simulation() {
    //std::cout << "Killing simulation\n";
    delete thePop;
    for (int si=0; si<sampleList.size(); ++si) {
        delete sampleList.at(si);
    }
}

void Simulation::runAndSave(int replicate) {
    thePop->initialize();
    // Sample generation 0:
    sampleList.push_back(thePop->makeSample());
    double t1 = getNow();
    double prevSec = 0;
    
    for(timeType t=1; t<=tmax; ++t) {
        //std::cout << "make gen " << t << '\n';
        thePop->makeNextGeneration();
        if (thePop->size()==0) {
            std::cout<< "Extinction at t = " << t << "!\n";
            break;
        }
        if( t % sampleInterval == 0 ) {
            sampleList.push_back(thePop->makeSample());
            if (verbose) {
                std::cout << "Generation " << t << ", population size " << thePop->size() << '\n';
                double nowSec = getNow()-t1;
                int timeLeft = int( (nowSec-prevSec)*(tmax-t)/sampleInterval );
                std::cout << "Time " << time2str(int(nowSec)) << ", time left : " << time2str(timeLeft) << "\n";
                prevSec = nowSec;
            }
        }
    }
    save(replicate);
}

void Simulation::save(int replicate) {
    int fileVersion = 1;
    std::ostringstream fileName;//(simName, std::ios_base::app);
    fileName << simName << "_results_" << replicate << ".sim";
    std::cout << "Saving " << fileName.str() << '\n';
    Simfile simfile(fileName.str());
    simfile.write(fileVersion);
    simfile.write((int)seed);
    simfile.writeString(parameterFileCopy.str());
    simfile.write(thePop->getGenetics().getLoci());
    // save all samples:
    simfile.write((int)sampleList.size());
    for(int si=0; si<sampleList.size(); ++si) {
        simfile << *sampleList.at(si);
    }
    if (thePop->isTrackingGenes()) {
        thePop->getGenetics().saveGeneLists(simfile);
    }
}
