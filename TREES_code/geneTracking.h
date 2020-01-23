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


class oSimfile;
class iSimfile;

class Gene {
    friend class GeneList;
    friend class Genetics;
protected:
    id_type id; // the id is assigned by the GeneList
    id_type parent;
    time_type birthTime, deathTime;
    float geneticEffect;
    std::vector<id_type> children;
    bool sampled; // sampled genes should not be pruned
    typedef std::vector<id_type>::iterator  childIter;
public:
    Gene(); // needed for lists, etc.
    Gene(id_type parent, time_type time, double effect);
    Gene(iSimfile& isf);
    void addChild(id_type child);
    void removeChild(id_type child);
    std::vector<id_type>& getChildren() { return children;}
    void writeToFile(oSimfile& osf);
};

//oSimfile& operator<<(oSimfile& os, Gene& a);

class GeneList : public std::vector<Gene> {
protected:
    id_type nextId;
public:
    GeneList();
    GeneList(iSimfile& isf);
    ~GeneList();
    id_type addGene(id_type parent, time_type time, double effect); // returns new Allele id
    Gene& getGene(id_type id);
    int getGeneIndex(id_type id);
    void pruneGenes(std::vector<bool>& alive);
    void writeGenes(oSimfile& osf);
};

//oSimfile& operator<<(oSimfile& os, GeneList& list);


#endif /* defined(__Species__geneTracking__) */
