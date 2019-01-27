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
    checkpointInterval = pf.getLong("checkpoint_interval"); // checkPoints==0 turns this feature off
    keepOldCheckpoints = pf.getBool("keep_old_checkpoints");
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

void Simulation::initialize(int replicate) {
    thePop->initialize();
    sampleList.clear();
    checkpointList.clear();
    makeFileName(replicate);
    // Sample generation 0:
    sampleList.push_back(thePop->makeSample());
}

void Simulation::makeFileName(int replicate) {
    std::ostringstream fileName;//(simName, std::ios_base::app);
    fileName << simName << "_results_" << replicate << ".sim";
    resultsFileName = fileName.str();
}

void Simulation::runFromGeneration(timeType gen) {
    double t1 = getNow();
    double prevSec = 0;
    
    for(timeType t=gen+1; t<=tmax; ++t) {
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
        if (checkpointInterval>0 && t%checkpointInterval==0) {
            unsigned seed = rand();
            if (keepOldCheckpoints) {
                checkpointList.push_back(thePop->makeCheckpoint(seed));
            } else {
                checkpointList.assign(1,thePop->makeCheckpoint(seed));
            }
            randseed(seed);
            save();
        }
    }
    if (checkpointInterval==0 || tmax%checkpointInterval != 0) {
        save();
    }
}


void Simulation::resume(std::string resultsFile, timeType resumeGen, std::string newResultsFile) {
    // read all samples and checkpoints:
    load(resultsFile);
    if (resumeGen==-1) {
        resumeGen = checkpointList.back()->getGeneration();
    } else {
        // Find the right checkpoint (they may not be saved at equal intervals):
        int cpi = 0;
        timeType gen = checkpointList[cpi]->getGeneration();
        while (gen<resumeGen && cpi<checkpointList.size()-1) {
            gen = checkpointList[++cpi]->getGeneration();
        }
        if (gen!=resumeGen) {
            std::cout << "Error. Can't find checkpoint for generation " << resumeGen << " in file " << resultsFile << '\n';
            exit(0);
        }
        // truncate lists:
        checkpointList.resize(cpi+1);
    }
    int si = 0;
    timeType gen = sampleList[si]->getGeneration();
    while (gen<resumeGen && si<sampleList.size()-1) {
        gen = sampleList[++si]->getGeneration();
    }
    if (gen>resumeGen) {
        --si;
    }
    sampleList.resize(si+1);
    thePop->resumeFromCheckpoint(*checkpointList.back());
    randseed(checkpointList.back()->seed);
    resultsFileName = newResultsFile;
    runFromGeneration(resumeGen);
}


void Simulation::save() {
    int fileVersion = 2; // As of v. 1.1
    std::cout << "Saving " << resultsFileName << '\n';
    oSimfile simfile(resultsFileName);
    simfile.write<int>(fileVersion);
    simfile.write<int>((int)seed);
    simfile.writeString(parameterFileCopy.str());
    simfile.write<int>(thePop->getGenetics().getLoci());
    // save all samples:
    simfile.write<int>((int)sampleList.size());
    for(Sample* s : sampleList) {
        thePop->writeSample(*s,simfile);
    }
    if (thePop->isTrackingGenes()) {
        thePop->getGenetics().writeGeneLists(simfile);
    }
    if (checkpointInterval>0) {
        simfile.write<int>((int)checkpointList.size());
        for(Checkpoint* cp : checkpointList) {
            thePop->writeCheckpoint(*cp, simfile);
        }
    }
}

// Load samples and checkpoints from results file. Don't change simulation parameters.
void Simulation::load(std::string& resultsFileName) {

    iSimfile resFile(resultsFileName);
    int fileVersion = resFile.read<int>();
    if( fileVersion!=2 ) {
        std::cout << "Wrong file version of results file " << resultsFileName << '\n';
        exit(0);
    }
    // read seed (not used):
    resFile.read<int>();
    // read (old) parameter-file copy (not used):
    resFile.readString();
    // read loci (not used):
    resFile.read<int>();
    
    int sampleCount = resFile.read<int>();
    sampleList.clear();
    sampleList.reserve(sampleCount);
    for (int si=0; si<sampleCount; ++si) {
        sampleList.push_back(thePop->readSample(resFile));
    }
    if (thePop->isTrackingGenes()) {
        thePop->getGenetics().readGeneLists(resFile);
    }

    if (checkpointInterval>0) {
        int count = resFile.read<int>();
        checkpointList.clear();
        checkpointList.reserve(count);
        for (int cpi=0; cpi<count; ++cpi) {
            checkpointList.push_back(thePop->readCheckpoint(resFile));
        }
    }

}

