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


#ifndef __Species__transform__
#define __Species__transform__

#include <stdio.h>

#include "parameterFile.h"
#include "types.h" // defines traitType
#include "fastmath.h"

class Genotype_Phenotype_map; // neede here to prevent circular header files.

class Transform {
public:
    virtual void transform(traitType* x, size_t count)=0;
    virtual double transform_single(double x)=0;
    virtual ~Transform();
};

class LinearTransform : public Transform {
protected:
    traitType offset;
    traitType scale;
public:
    LinearTransform(ParameterFile& pf);
    LinearTransform(traitType offset, traitType scale);
    virtual ~LinearTransform();
    virtual void transform(traitType* x, size_t count);
    virtual double transform_single(double x);
};


class Range : public LinearTransform {
public:
    Range(ParameterFile& pf, Genotype_Phenotype_map& GP_map);
    virtual ~Range();
};


class AbsTransform : public Transform {
public:
    AbsTransform(ParameterFile& pf);
    virtual ~AbsTransform();
    virtual void transform(traitType* x, size_t count);
    virtual double transform_single(double x);
};

class LogisticTransform : public Transform {
    traitType xmin;
    traitType scale;
    Fastexp fexp;
public:
    LogisticTransform(ParameterFile& pf);;
    virtual ~LogisticTransform();
    virtual void transform(traitType* x, size_t count);
    virtual double transform_single(double x);
};

class NormalDeviate : public Transform {
    traitType SD;
public:
    NormalDeviate(ParameterFile& pf);
    virtual ~NormalDeviate();
    virtual void transform(traitType* x, size_t count);
    virtual double transform_single(double x);
};

#endif /* defined(__Species__transform__) */
