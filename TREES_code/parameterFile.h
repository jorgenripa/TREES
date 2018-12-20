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


#ifndef __Species__parameterFile__
#define __Species__parameterFile__

#include <stdio.h>
#include <cstdio>
#include <string>
#include <fstream>
#include <cstdint>

class ParameterFile {
public:
    ParameterFile(const std::string& fileName);
    ~ParameterFile();
    double getDouble(const char* name);
    double getPositiveDouble(const char* name);
    int64_t getLong(const char* name);
    int64_t getPositiveLong(const char* name);
    int getInt(const char* name) { return (int)getLong(name); }
    int getPositiveInt(const char* name);
    bool getBool(const char* name);
    std::string getString(const char* name);
    std::string getStringLower(const char* name);
    void getStringPair(std::string& s1, std::string& s2);
    void close();
    bool eof() { return pfile.eof(); }
protected:
    std::ifstream pfile;
    char skipWhiteSpace();
    void skipRestOfLine();
    void getNextLine(char* buffer);
};


#endif /* defined(__Species__parameterFile__) */
