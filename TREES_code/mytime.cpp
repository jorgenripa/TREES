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


#include <ctime>
#include <sstream>
#include <cmath>

#include "mytime.h"


double getNow() {
    return clock()/double(CLOCKS_PER_SEC);
}

std::string time2str(double secs) {
    std::stringstream s;
    s << int(secs/3600) << "h " << int(secs/60) % 60 << "m " << fmod(secs,60.0) << 's';
    return s.str();
}
