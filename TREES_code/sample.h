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


#ifndef __Species__sample__
#define __Species__sample__

#include <stdio.h>
#include <vector>

#include "simfile.h"
#include "types.h"

class SampleData {
    // Contains one data point per individual
protected:
    SampleData();
public:
    virtual ~SampleData();
    virtual void writeToFile(Simfile& os) = 0;
};

class FloatData : public SampleData {
protected:
    std::vector<float> theData;
public:
    FloatData(std::vector<double>& data);
    virtual ~FloatData();
    virtual void writeToFile(Simfile& os);
};

class IntData : public SampleData {
protected:
    std::vector<int> theData;
public:
    virtual ~IntData();
    IntData(std::vector<int>& data);
    virtual void writeToFile(Simfile& os);
};

class IdData : public SampleData {
protected:
    std::vector<idType> theData;
public:
    virtual ~IdData();
    IdData(std::vector<idType>& data);
    virtual void writeToFile(Simfile& os);
};

class Sample {
    friend class Population;
    friend  Simfile& operator<<(Simfile& os, Sample& sample);
protected:
    timeType generation;
    int individualCount;
    std::vector<SampleData*> dataList;
public:
    Sample(timeType generation, int n);
    ~Sample();
    void addData(SampleData* sd);
};

Simfile& operator<<(Simfile& os, Sample& sample);


#endif /* defined(__Species__sample__) */
