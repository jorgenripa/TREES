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


#include <cstdlib>
#include <cctype>
#include <cstring>
#include <iostream>
#include <limits>
#include <algorithm>

#include "parameterFile.h"

ParameterFile::ParameterFile(const std::string& fileName) :
pfile(fileName)
{
    if (!pfile.is_open()) {
        std::cout << "Can't open file " << fileName << "!\n";
        exit(0);
    }
}

ParameterFile::~ParameterFile() {
    //std::cout << "Killing pf\n";
    close();
}

void ParameterFile::close() {
    if (pfile.is_open()) {
        pfile.close();
    }
}

// a few ad-hoc low-level functions:
bool my_isspace(const char c) {
    return isspace(c) || c=='\xa0'; // add non-breaking space to the list
}

int skipWhiteSpace(const char* str, int pos) {
    while (str[pos]!=0 && my_isspace(str[pos])) {
        ++pos;
    }
    return pos;
}

int getString(const char* instr, int pos, std::string& S) {
    char buffer[100];
    buffer[0] = 0;
    int bufpos = 0;
    while ( instr[pos]>0 && !my_isspace(instr[pos]) && !(instr[pos]==':') && bufpos<99) {
        buffer[bufpos++] = instr[pos++];
    }
    buffer[bufpos] = 0;
    if (bufpos>=99) {
        std::cout << "Buffer overflow : " << buffer << '\n';
        exit(0);
    }
    S = buffer;
    return pos;
}

void stringToLower( std::string& s ) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
}

void ParameterFile::getStringPair(std::string& s1, std::string& s2) {
    char line[500];
    getNextLine(line);
    int pos = 0;
    pos = ::skipWhiteSpace(line,pos);
    if (line[pos]==0) {
        s1 = "";
        s2 = "";
    } else {
        pos = ::getString(line,pos,s1);
        pos = ::skipWhiteSpace(line,pos);
        if (line[pos]==':') {
            ++pos;
        } else {
            std::cout << "Bad input : " << line << '\n';
            exit(0);
        }
        pos = ::skipWhiteSpace(line,pos);
        pos = ::getString(line,pos,s2);
    }
    stringToLower(s1);
    //stringToLower(s2);
}

double ParameterFile::getDouble(const char* name) {
    std::string s1, s2;
    getStringPair(s1,s2);
    if (s1 != name ) {
        std::cout << "Parameter error : \n";
        std::cout << "Reading parameter : " << name << '\n';
        std::cout << "Found parameter : " << s1 << '\n';
        exit(1);
    }
    return std::atof(s2.c_str());
}

int64_t ParameterFile::getLong(const char* name) {
    // read a double instead of digit to interpret input as "1e6" correctly:
    return (int64_t)getDouble(name);
}

int64_t ParameterFile::getPositiveLong(const char* name) {
    int64_t i64 = getLong(name);
    if (i64<=0) {
        std::cout << "Parameter " << name << " has to be positive. Found value : " << i64 << '\n';
        exit(0);
    }
    return i64;
}

double ParameterFile::getPositiveDouble( const char* name) {
    double d = getDouble(name);
    if (d<=0) {
        std::cout << "Parameter " << name << " has to be positive. Found value : " << d << '\n';
        exit(0);
    }
    return d;
    
}

int ParameterFile::getPositiveInt(const char* name) {
    int i = (int)getLong(name);
    if (i<=0) {
        std::cout << "Parameter " << name << " has to be positive. Found value : " << i << '\n';
        exit(0);
    }
    return i;
}

bool ParameterFile::getBool(const char* name) {
    std::string boolstr = getString(name);
    bool value;
    if (std::toupper(boolstr[0])=='Y') {
        value = true;
    } else if (std::toupper(boolstr[0])=='N') {
        value = false;
    } else {
        std::cout << "Y/N parameter error : \n";
        std::cout << "Reading parameter : " << name << '\n';
        std::cout << "Found value : " << boolstr << '\n';
        exit(1);
    }
    return value;
}


std::string ParameterFile::getString(const char* name) {
    std::string s1, s2;
    getStringPair(s1,s2);
    if (s1!=name ) {
        std::cout << "Parameter error : \n";
        std::cout << "Reading parameter : " << name << '\n';
        std::cout << "Found parameter : " << s1 << '\n';
        exit(1);
    }
    return s2;
}

std::string ParameterFile::getStringLower(const char* name) {
    std::string s(getString(name));
    stringToLower(s);
    return s;
}

char ParameterFile::skipWhiteSpace() {
    char c = pfile.get();
    while (my_isspace(c)) { // skips space, form feed, line feed, carriage return, tab + nonbreaking space
        c = pfile.get();
    }
    return c;
}

void ParameterFile::skipRestOfLine() {
    char lastc = pfile.get();
    while (!(pfile.eof() || lastc=='\n' || lastc=='\r' )) {
        lastc = pfile.get();
    }
}

void ParameterFile::getNextLine(char *buffer) {
    int lastpos;
    // Skip comment lines starting with #, and blank lines
    buffer[0] = skipWhiteSpace();
    while (buffer[0]=='#') {
        skipRestOfLine();
        buffer[0] = skipWhiteSpace();
    }
    if (pfile.eof()) {
        buffer[0] = 0;
        return;
    }
    lastpos = 0;
    // read until comma, #, linefeed, or EOF:
    while (!(pfile.eof() || buffer[lastpos]=='\n' || buffer[lastpos]=='\r' || buffer[lastpos]==',' || buffer[lastpos]=='#'  )) {
        char lastc = pfile.get();
        buffer[++lastpos] = lastc;
    }
    // skip comment at end of line:
    if (buffer[lastpos]=='#') {
        skipRestOfLine();
    }
    buffer[lastpos] = 0;
}
