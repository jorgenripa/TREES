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
    
    // Check command line parameters: the parameter file and possible replicate indeces
    std::string parameterFileName;
    int rep1, rep2;
    if (argc<2) {
        std::cout << "Error: Missing parameter file. \n Syntax:\n Species parFile.txt [rep1] [rep2]\n";
        return 0;
    } else {
        parameterFileName = argv[1];
        if(argc>=3) {
            rep1 = atoi(argv[2]);
            if(argc>=4) {
                rep2 = atoi(argv[3]);
            } else {
                rep2 = rep1;
            }
        } else {
            rep1 = 1;
            rep2 = 1;
        }
    }

    // main loop running replicates:
    std::cout << "Running replicate " << rep1 << " to " << rep2 << '\n';
    for (int i=rep1; i<=rep2; ++i) {
        Simulation theSim(parameterFileName);
        double t1 = getNow();
        std::cout << "Running " << i << '\n';
        theSim.runAndSave(i);
        std::cout << time2str(getNow()-t1) << '\n';
    }
    return 0;    
}
