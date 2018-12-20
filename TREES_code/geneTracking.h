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


#ifndef __Species__geneTracking__
#define __Species__geneTracking__

#include <stdio.h>
#include <vector>
#include <cstdint>
#include "types.h"

typedef uint64_t idType;

class Simfile;

class Gene {
    friend class GeneList;
    friend class Genetics;
protected:
    idType id; // the id is assigned by the GeneList
    idType parent;
    timeType birthTime, deathTime;
    double geneticEffect;
    std::vector<idType> children;
    bool sampled; // sampled genes should not be pruned
    typedef std::vector<idType>::iterator  childIter;
public:
    Gene(); // somehow, this is needed
    Gene(idType parent, timeType time, double effect);
    void addChild(idType child);
    void removeChild(idType child);
    std::vector<idType>& getChildren() { return children;}
    friend  Simfile& operator<<(Simfile& os, Gene& a);
};

Simfile& operator<<(Simfile& os, Gene& a);

class GeneList : public std::vector<Gene> {
protected:
    idType nextId;
public:
    GeneList();
    ~GeneList();
    idType addGene(idType parent, timeType time, double effect); // returns new Allele id
    Gene& getGene(idType id);
    int getGeneIndex(idType id);
    void pruneGenes(std::vector<bool>& alive);
};

#endif /* defined(__Species__geneTracking__) */
