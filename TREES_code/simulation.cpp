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
#include "mytime.h"

Simulation::Simulation(std::string& parameterFileName, bool verbose){
    ParameterFile pf(parameterFileName);

    this->verbose = verbose;
    simName = parameterFileName;
    // skip file extension:
    simName.erase(simName.find_last_of('.'));
    std::cout << "name: " << simName << '\n';
    tmax = pf.getPositiveLong("t_max");
    sampleInterval = pf.getPositiveLong("sample_interval");
    std::string option_string = pf.getStringLower("microsamples");
    microsample_option = option_string[0];
    checkpointInterval = pf.getLong("checkpoint_interval"); // checkPoints==0 turns this feature off
    keepOldCheckpoints = pf.getBool("keep_old_checkpoints");
    std::string seedS = pf.getStringLower("seed");
    if (seedS[0]=='r') {
        seed = randomize();
    } else {
        seed = (seed_type)std::atoll(seedS.c_str());
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
    for (size_type si=0; si<sample_list.size(); ++si) {
        delete sample_list.at(si);
    }
    for (size_type si=0; si<microsample_list.size(); ++si) {
        delete microsample_list.at(si);
    }
}

void Simulation::initialize(int replicate) {
    thePop->initialize();
    sample_list.clear();
    microsample_list.clear();
    checkpointList.clear();
    makeFileName(replicate);
    // Sample generation 0:
    sample_list.push_back(new Population_sample(*thePop,0));
    if (microsample_option!='n') {
        microsample_list.push_back(new Microsample(*thePop,microsample_option,0));
    }
}

void Simulation::makeFileName(int replicate) {
    std::ostringstream fileName;//(simName, std::ios_base::app);
    fileName << simName << "_results_" << replicate << ".sim";
    resultsFileName = fileName.str();
}

void Simulation::runFromGeneration(time_type gen, seed_type gen_seed, double time_so_far) {
    double start_time = getNow();
    double prevSec = time_so_far;
    randseed(gen_seed);
    for(time_type t=gen+1; t<=tmax; ++t) {
        if (verbose) {
            std::cout << t << std::flush;
        }
        //std::cout << "make gen " << t << '\n';
        thePop->makeNextGeneration();
        // Seed for next generation (may be saved in checkpoint):
        seed_type next_seed = make_seed();
        randseed(next_seed);

        double cputime = getNow()-start_time+time_so_far;
        
        if (verbose) {
            int backsteps = int(log10(t))+1;
            std::cout << std::string(backsteps,'\b');
        }
        if (thePop->size()==0) {
            std::cout<< "Extinction at t = " << t << "!\n";
            break;
        }
        if (microsample_option!='n') {
            microsample_list.push_back(new Microsample(*thePop, microsample_option, cputime));
        }
        if( t % sampleInterval == 0 ) {
            sample_list.push_back(new Population_sample(*thePop,cputime));
            if (verbose) {
                std::cout << "Generation " << t << ", population size " << thePop->size() << '\n';
                int timeLeft = int( (cputime-prevSec)*(tmax-t)/sampleInterval );
                std::cout << "Time " << time2str(int(cputime)) << ", time left : " << time2str(timeLeft) << "\n";
                prevSec = cputime;
            }
        }
        if (checkpointInterval>0 && t%checkpointInterval==0) {
            Population_checkpoint* cp = new Population_checkpoint(*thePop,next_seed,cputime);
            if (keepOldCheckpoints) {
                checkpointList.push_back(cp);
            } else {
                checkpointList.assign(1,cp);
            }
            save();
        }
    }
    if (checkpointInterval==0 || tmax%checkpointInterval != 0) {
        save();
    }
}

void truncate_sample_list(std::vector<Sample_base*>& list, time_type last_gen) {
    if (list.size()>0) {
        int si=0;
        time_type gen = list.at(si)->get_generation();
        while (gen<last_gen && si<list.size()-1) {
            gen = list.at(++si)->get_generation();
        }
        if (gen>last_gen) {
            --si;
        }
        // si is now index of last generation, si+1 is size of valid list
        if (list.size() > si+1) {
            for (int del_i=si+1; del_i<list.size(); ++del_i) {
                delete list.at(del_i);
            }
            list.resize(si+1);
        }
    }
}

void Simulation::resume(std::string resultsFile, time_type resumeGen, std::string newResultsFile) {
    // read all samples and checkpoints:
    load(resultsFile);
    if (resumeGen==-1) { // use last checkpoint:
        resumeGen = checkpointList.back()->get_generation();
    } else {
        // Find the right checkpoint (they may not be saved at equal intervals):
        truncate_sample_list(checkpointList, resumeGen);
        if (checkpointList.size()==0 || checkpointList.back()->get_generation() !=resumeGen) {
            std::cout << "Error. Can't find checkpoint for generation " << resumeGen << " in file " << resultsFile << '\n';
            exit(0);
        }
    }
    // truncate microsample list:
    truncate_sample_list(microsample_list, resumeGen);
    // truncate sample list:
    truncate_sample_list(sample_list, resumeGen);

    Population_checkpoint& cp = *dynamic_cast<Population_checkpoint*>(checkpointList.back());
    thePop->resumeAtCheckpoint(cp);
    resultsFileName = newResultsFile;
    runFromGeneration(resumeGen, cp.get_seed(), cp.get_cputime());
}


void Simulation::save() {
    std::cout << "Saving " << resultsFileName << "..." << std::flush;
    oSimfile simfile(resultsFileName);
    simfile.write<seed_type>(seed);
    simfile.writeString(parameterFileCopy.str());
    simfile.write<size_type>(thePop->getGenetics().getLoci());
    // save microsamples
    simfile.write<size_type>((size_type)microsample_list.size());
    for(Sample_base* ms : microsample_list) {
        dynamic_cast<Microsample*>(ms)->write_to_file(simfile);
    }
    // save samples:
    simfile.write<size_type>((size_type)sample_list.size());
    for(Sample_base* s : sample_list) {
        dynamic_cast<Population_sample*>(s)->write_to_file(simfile);
    }
    if (thePop->isTrackingGenes()) {
        thePop->getGenetics().writeGeneLists(simfile);
    }
    if (checkpointInterval>0) {
        simfile.write<size_type>((size_type)checkpointList.size());
        for(Sample_base* cp : checkpointList) {
            dynamic_cast<Population_checkpoint*>(cp)->write_to_file(simfile);
        }
    }
    simfile.close();
    std::cout << "Done.\n" << std::flush;;
}

// Load samples and checkpoints from results file. Don't change simulation parameters.
void Simulation::load(std::string& resultsFileName) {

    iSimfile resFile(resultsFileName);
    // read seed (this should be kept):
    seed = resFile.read<seed_type>();
    // read (old) parameter-file copy (not used):
    resFile.readString();
    // read loci (not used):
    resFile.read<size_type>();
    
    // read microsamples
    size_type microsample_count = resFile.read<size_type>();
    microsample_list.clear();
    for (size_type si=0; si<microsample_count; ++si) {
        microsample_list.push_back(new Microsample(resFile));
    }
    // Read samples:
    size_type sampleCount = resFile.read<size_type>();
    sample_list.clear();
    sample_list.reserve(sampleCount);
    for (size_type si=0; si<sampleCount; ++si) {
        sample_list.push_back(new Population_sample(resFile));
    }
    if (thePop->isTrackingGenes()) {
        thePop->getGenetics().readGeneLists(resFile);
    }

    if (checkpointInterval>0) {
        size_type count = resFile.read<size_type>();
        checkpointList.clear();
        checkpointList.reserve(count);
        for (size_type cpi=0; cpi<count; ++cpi) {
            checkpointList.push_back(new Population_checkpoint(resFile));
        }
    }

}

