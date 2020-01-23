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


#include <cmath>
#include <iostream>

#include "transform.h"
#include "genotype_phenotype_map.h"
#include "fastmath.h"
#include "parameterFile.h"
#include "randgen.h"
#include "trait.h"
#include "genetics.h"

Transform::~Transform() {}

LinearTransform::LinearTransform(ParameterFile& pf) {
    offset = pf.getDouble("offset");
    scale = pf.getDouble("scale");
}

LinearTransform::LinearTransform(traitType offset, traitType scale) {
    this->offset = offset;
    this->scale = scale;
}

LinearTransform::~LinearTransform() {}

void LinearTransform::transform(traitType* x, size_t count) {
    for (size_t i=0; i<count; ++i) {
        *x = offset + scale*(*x);
        ++x;
    }
}

double LinearTransform::transform_single(double x) {
    return offset + scale*x;
}

Range::Range(ParameterFile& pf, Genotype_Phenotype_map& GP_map): LinearTransform(0,0) {
    Diallelic_GP_map* the_map = dynamic_cast<Diallelic_GP_map*>(&GP_map);
    if (the_map) {
        traitType min_value = pf.getDouble("min");
        traitType max_value = pf.getDouble("max");
        int L = the_map->get_loci_per_dim(); // untransformed values [-L,L]
        offset = L + min_value;
        scale = (max_value-min_value)/2.0/L;
    } else {
        std::cout << "Error! \n" << "Trait " << GP_map.get_trait().getName() << ": A Range transform requires Diallelic Genetics.\n";
        exit(1);
    }
}

Range::~Range() {};

AbsTransform::AbsTransform(ParameterFile& pf) {}
AbsTransform::~AbsTransform() {}

void AbsTransform::transform(traitType* x, size_t count) {
    for (size_t i=0; i<count; ++i) {
        *x = fabs(*x);
        ++x;
    }
}
double AbsTransform::transform_single(double x) {
    return fabs(x);
}


LogisticTransform::LogisticTransform(ParameterFile& pf) {
    xmin = pf.getDouble("min");
    double xmax = pf.getDouble("max");
    scale = xmax - xmin;
}

LogisticTransform::~LogisticTransform() {}

void LogisticTransform::transform(traitType* x, size_t count){
    // A logistic transform to the range (xmin, xmax)
    // x=0 is mapped to the range midpoint
    // The slope at x=0 is scaled to 1
    // f = xmin + scale/(1 + exp(-4x/scale))
    // where scale = xmax-xmin
    for (size_t i=0; i<count; ++i) {
        double expmx = fexp.exp(-4*(*x)/scale);
        *x = xmin + scale/(1+expmx);
        ++x;
    }
}

double LogisticTransform::transform_single(double x) {
    double expmx = fexp.exp(-4*x/scale);
    x = xmin + scale/(1+expmx);
    return x;
}


NormalDeviate::NormalDeviate(ParameterFile& pf) {
    SD = pf.getPositiveDouble("sd");
}

NormalDeviate::~NormalDeviate() {}

void NormalDeviate::transform(traitType* x, size_t count) {
    for (size_t i=0; i<count; ++i) {
        *x += SD*randn();
        ++x;
    }
}
double NormalDeviate::transform_single(double x) {
    return x + SD*randn();
}
