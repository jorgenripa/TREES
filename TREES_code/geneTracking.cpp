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
#include "stdlib.h"

#include "geneTracking.h"
#include "simfile.h"
#include "types.h"

Gene::Gene(idType parent, timeType time, double effect) {
    id = 0;
    this->parent = parent;
    birthTime = time;
    deathTime = -1;
    geneticEffect = effect;
    children.reserve(10);
    sampled = false;
}

Gene::Gene() {
    id = 0;
    parent = 0;
    birthTime = 0;
    deathTime = -1;
    geneticEffect = 0;
    children.reserve(10);
    sampled = false;
}

void Gene::addChild(idType child) {
    children.push_back(child);
}

void Gene::removeChild(idType child) {
    for (childIter c=children.begin(); c!=children.end(); ++c) {
        if (*c==child) {
            children.erase(c);
            return;
        }
    }
    std::cout << "removeChild error, child not found!\n";
    exit(1);
}

Simfile& operator<<(Simfile& sf, Gene& g) {
    sf.write(g.id);
    sf.write(g.parent);
    sf.write(g.birthTime);
    sf.write(g.deathTime);
    sf.write((float)g.geneticEffect);
    sf.write((int)g.children.size());
    if (g.children.size()>0) {
        sf.writeArray(&g.children[0], (int)g.children.size());
    }
    return sf;
}


GeneList::GeneList() {
    clear();
    reserve(10000);
    nextId = 1; // 0 is reserved for invalid Genes
}

GeneList::~GeneList() {
}

idType GeneList::addGene(idType parent, timeType time, double effect) {
    Gene newA( parent, time, effect);
    newA.id = nextId++;
    push_back(newA);
    if(parent>0) {
        getGene(newA.parent).addChild(newA.id);
    }
    return newA.id;
}

int GeneList::getGeneIndex(idType id) {
    for (int a=0; a<size(); ++a) {
        if (at(a).id==id) {
            return a;
        }
    }
    return -1;
    // Comment: This function does not return an error, since it is sometimes used to test
    // the presence of a particular id in the list.
    // Care need to be taken when using the result.
}

Gene& GeneList::getGene(idType id) {
    int i = getGeneIndex(id);
    if (i>=0) {
        return at(getGeneIndex(id));
    } else {
        std::cout << "getGene error, id not found!\n";
        exit(1);
    }
}

// remove and delete all non-sampled Genes:
void GeneList::pruneGenes(std::vector<bool>& alive) {
    int iw = 0;
    for (int ir=0; ir<size(); ++ir) {
        if(at(ir).sampled || alive.at(ir)) {
            if (ir>iw) {
                at(iw) = at(ir);
            }
            ++iw;
        } else { // remove dead gene from parent's list of children
            idType parent = at(ir).parent;
            if (parent>0) { // zero is the parent of the root
                int parent_i = getGeneIndex(parent);
                // If parent still in the list and active:
                if(parent_i>=0 && (at(parent_i).sampled || alive.at(parent_i))) {
                    at(parent_i).removeChild(at(ir).id);
                }
            }
        }
    }
    resize(iw);
}
