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



#ifndef __Species__fastexp__
#define __Species__fastexp__

#include <stdio.h>
#include <cmath>


class Fastexp {
    double* etable;
    double* edifftable;
public:
    static const int MAXX;
    Fastexp();
    // This routine has a max relative error of approx 0.002
    // Relies on fillExpTable() above
    inline double exp(double x) {
        if (x==0) {
            return 1.0;
        } else if (x>0) {
            // avoid recursive inline
            if (x>MAXX) {
                return std::exp(x);
            }
            double z = x*10;
            int i = (int)z;
            double y1 = etable[i];
            return 1.0/(y1 + (z-i)*edifftable[i]);
        } else if (x<-MAXX) {
            return std::exp(x);
        } else {
            double z = -x*10;
            int i = (int)z;
            double y1 = etable[i];
            return y1 + (z-i)*edifftable[i];
        }
    }
};


// This HAS TO be called once before calling fastexp2:
void fillExpTable();

inline double absipow(double base, int exp)
{
    if (exp<0) {
        return 1.0/absipow(base,-exp);
    }
    if (base<0)
        base = -base;
    switch (exp) {
        case 0:
            return 1;
        case 1:
            return base;
        case 2:
            return base*base;
        default:
            int result = 1;
            while (exp)
            {
                if (exp & 1)
                    result *= base;
                exp >>= 1;
                base *= base;
            }
            
            return result;
    }
}

// This is used by Fitness module StabilizingSelection:
inline double fastabspow(double base, double exp) {
    if (base<0) {
        base = -base;
    }
    if (exp-(int)exp == 0) {
        return absipow(base,(int)exp);
    } else {
        return pow(base,exp);
    }
}

#endif /* defined(__Species__fastexp__) */
