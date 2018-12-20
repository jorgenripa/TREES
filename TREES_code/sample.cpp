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


#include <iostream>

#include "sample.h"
#include "simfile.h"

SampleData::SampleData()
{}

SampleData::~SampleData()
{}

FloatData::FloatData(std::vector<double>& data) {
    // copy data to theData, converting to float
    theData.resize(data.size());
    for (int i=0; i<data.size(); ++i) {
        theData.at(i) = (float)data.at(i);
    }
}

FloatData::~FloatData()
{}

void FloatData::writeToFile(Simfile& os) {
    os.writeArray(&theData.at(0), (int)theData.size());
}

IntData::~IntData()
{}

IntData::IntData(std::vector<int>& data) {
    theData = data;
}

void IntData::writeToFile(Simfile& os) {
    os.writeArray(&theData.at(0), (int)theData.size());
}

IdData::IdData(std::vector<idType>& data) {
    theData = data;
}

IdData::~IdData()
{}

void IdData::writeToFile(Simfile& os) {
    os.writeArray(&theData.at(0), (int)theData.size());
}


Sample::Sample(timeType generation, int n) {
    this->generation = generation;
    individualCount = n;
    dataList.clear();
}

Sample::~Sample() {
    for( SampleData* sdp : dataList) {
        delete sdp;
    }
}

void Sample::addData(SampleData* sd) {
    dataList.push_back(sd);
}

Simfile& operator<<(Simfile& os, Sample& sample) {
    os.write(sample.generation);
    os.write(sample.individualCount);
    for( SampleData* sdp : sample.dataList) {
        sdp->writeToFile(os);
    }
    return os;
}
