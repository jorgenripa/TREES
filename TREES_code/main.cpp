//
//  main.cpp
//  TREES project
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
#include <iostream>
#include <string>
#include <locale>

#include "fastmath.h"
#include "parameterFile.h"
#include "simulation.h"
#include "mytime.h"

// Main procedure
int main(int argc, const char * argv[]) {

    // Calculate table for interpolation of exponential function:
    fillExpTable();// don't remove!!!!
    
    int next_par_i = 1;
    // check for command line options:
    bool resume = false;
    bool verbose = false;
    while (next_par_i<argc && argv[next_par_i][0]=='-') {
        std::string option(argv[next_par_i]);
        if (option=="-resume") {
            resume = true;
        } else if (option=="-v") {
            verbose = true;
        } else {
            std::cout << "Unknown option : " << argv[next_par_i] << '\n';
            exit(1);
        }
        ++next_par_i;
    }
    // Check command line parameters: the parameter file and possible replicate indeces
    // syntax: > TREES parameterFile.txt [rep1 [rep2]]
    std::string parameterFileName;
    int rep1, rep2;
    if (argc-next_par_i < 1) {
        std::cout << "Error: Missing parameter file. \n Syntax:\n TREES parameter_file [first_replicate  [last_replicate]]\n";
        return 1;
    }
    parameterFileName = argv[next_par_i++];
    if (resume) {
        // resume option, single replicate
        // syntax: > TREES -resume parameter_file.txt results_file.txt [generation [new_results_file.txt]]
        std::string resultsFile, newResultsFile;
        time_type gen;
        if (next_par_i<argc) {
            resultsFile = argv[next_par_i++];
            if (next_par_i<argc) {
                // read a float instead of digit to interpret input as "1e6" correctly:
                gen = atof(argv[next_par_i++]);
            } else {
                gen = -1; // meaning: start at last available checkpoint
            }
            Simulation theSim(parameterFileName, verbose);
            if (next_par_i<argc) {
                newResultsFile = argv[next_par_i++];
            } else {
                newResultsFile = resultsFile;
            }
            theSim.resume(resultsFile,gen,newResultsFile);
        } else {
            std::cout << "Error: Missing resume parameters.\n Syntax:\n TREES [-resume] [-v] parameter_file results_file [generation [new_results_file]]\n Omitting the generation means starting at last saved checkpoint.\n" ;
            return 1;
        }
    } else {
        // standard run option, one or more replicates
        if(next_par_i<argc) {
            rep1 = atoi(argv[next_par_i++]);
            if(next_par_i<argc) {
                rep2 = atoi(argv[next_par_i]);
            } else {
                rep2 = rep1;
            }
        } else {
            rep1 = 1;
            rep2 = 1;
        }
        
        
        // main loop running replicates:
        std::cout << "Running replicate " << rep1 << " to " << rep2 << '\n';
        Simulation theSim(parameterFileName, verbose);
        for (int i=rep1; i<=rep2; ++i) {
            std::cout << "Running " << i << '\n';
            theSim.initialize(i);
            double t1 = getNow();
            theSim.runFromGeneration(0, theSim.get_seed(), 0);
            std::cout << time2str(getNow()-t1) << '\n';
        }
    }
    return 0;
}
