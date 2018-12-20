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
#include "transform.h"
#include "fastmath.h"
#include "parameterFile.h"
#include "randgen.h"

Transform::~Transform() {}

LinearTransform::LinearTransform(ParameterFile& pf) {
    offset = pf.getDouble("offset");
    scale = pf.getDouble("scale");
}
LinearTransform::~LinearTransform() {}

void LinearTransform::transform(double* x, size_t count) {
    for (size_t i=0; i<count; ++i) {
        *x = offset + scale*(*x);
        ++x;
    }
}

AbsTransform::AbsTransform(ParameterFile& pf) {}
AbsTransform::~AbsTransform() {}

void AbsTransform::transform(double* x, size_t count) {
    for (size_t i=0; i<count; ++i) {
        *x = fabs(*x);
        ++x;
    }
}

LogisticTransform::LogisticTransform(ParameterFile& pf) {
    xmin = pf.getDouble("min");
    xmax = pf.getDouble("max");
}

LogisticTransform::~LogisticTransform() {}

void LogisticTransform::transform(double* x, size_t count){
    // A logistic transform to the range (xmin, xmax)
    // x=0 is mapped to the range midpoint
    // The slope at x=0 is scaled to 1
    // f = xmin + scale/(1 + exp(-4x/scale))
    // where scale = xmax-xmin
    Fastexp fexp;
    double scale = xmax-xmin;
    for (size_t i=0; i<count; ++i) {
        double expmx = fexp.exp(-4*(*x)/scale);
        *x = xmin + scale/(1+expmx);
        ++x;
    }
}

NormalDeviate::NormalDeviate(ParameterFile& pf) {
    SD = pf.getPositiveDouble("sd");
}

NormalDeviate::~NormalDeviate() {}

void NormalDeviate::transform(double* x, size_t count) {
    for (size_t i=0; i<count; ++i) {
        *x += SD*randn();
        ++x;
    }
}
