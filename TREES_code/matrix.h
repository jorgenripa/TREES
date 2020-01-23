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

#ifndef matrix_h
#define matrix_h

#include <vector>
#include <cstring>
#include "simfile.h"

// Matrix class without default constructing
template<typename T> class Matrix {
protected:
    T* X;
    uint32_t M, N, X_size;
    void initialize() {M=0; N=0; X=NULL; X_size=0;}
public:
    Matrix() {initialize();}

    void reserve( uint32_t Mres, uint32_t Nres) {
        if (Mres*Nres > X_size) {
            T* newX = new T[Mres*Nres];
            if (X && M*N>0) {
                std::memcpy(newX, X, M*N*sizeof(T));
            }
            if (X) {
                delete [] X;
            }
            X = newX;
            X_size = Mres*Nres;
        }
    }
    
    void assign(uint32_t newM, uint32_t newN, T val) {
        reserve(newM,newN);
        M = newM;
        N = newN;
        T* i = X;
        for (int c=0; c<N; ++c) {
            for (int r=0; r<M; ++r) {
                *i= val;
                ++i;
            }
        }
    }
    
    void resize(uint32_t newM, uint32_t newN) {
        reserve(newM,newN);
        M=newM; N = newN;
    }

    Matrix(uint32_t M, uint32_t N) { initialize(); resize(M, N);}
    Matrix(uint32_t M, uint32_t N, T val) { initialize(); assign(M,N,val); }
    // copy constructor:
    Matrix(const Matrix<T>& M2) {
        initialize();
        *this = M2;
    }
    
    // assignment operator:
    Matrix<T>& operator = (const Matrix<T> &M2)
    {
        resize(M2.M, M2.N);
        std::memcpy(X, M2.X, M2.M*M2.N*sizeof(T));
        return *this;
    }
    
    Matrix(iSimfile& isf) {
        initialize();
        M = isf.read<uint32_t>();
        N = isf.read<uint32_t>();
        reserve(M,N);
        isf.readCArray<T>(X, M*N);
    }
    ~Matrix() { if (X) delete [] X; }
    
    void write_to_file(oSimfile& osf) {
        osf.write<uint32_t>(M);
        osf.write<uint32_t>(N);
        osf.writeCArray<T>(X, M*N);
    }
    void read_from_file(iSimfile & isf) {
        M = isf.read<uint32_t>();
        N = isf.read<uint32_t>();
        reserve(M,N);
        isf.readCArray<T>(X, M*N);
    }
    T& operator() (uint32_t row, uint32_t col) {
        if (row<M && col<N) {
            return X[row+col*M];
        }
        else {
            std::cout << "Index out of bounds : " << row << ", " << col << '\n';
            std::cout << "M = " << M << ", N = " << N << '\n';
            exit(9);
        }
    }
    uint32_t get_M() { return M;}
    uint32_t get_N() { return N;}
    const T* getX() { return X;}
    uint32_t capacity() { return X_size; }
    
    void add_column(T* values) {
        resize(M,N+1);
        T* new_pos = &(*this)(0,N-1);
        for (uint32_t d=0; d<M; ++d) {
            *new_pos = *values;
            ++new_pos;
            ++values;
        }
    }
    double row_mean(uint32_t row){
        double sum=0;
        T* xp = &(*this)(row,0);
        for (uint32_t col=0; col<N; ++col) {
            sum += *xp;
            xp += M;
        }
        return sum/N;
    }
    double col_mean(uint32_t col){
        double sum=0;
        T* xp = &(*this)(0,col);
        for (uint32_t row=0; row<M; ++row) {
            sum += *xp;
            ++xp;
        }
        return sum/M;
    }
    double row_variance(uint32_t row, double row_mean) {
        // Use the row_mean for better precision than calculating the total sum of squares separately
        double sqsum = 0;
        T* xp = &(*this)(row,0);
        for (int col=0; col < N; ++col) {
            double dx = *xp - row_mean; // better precision than
            sqsum += dx*dx;
            xp += M;
        }
        return sqsum/N;
    }

    void compact_data(std::vector<bool>& alive) {
        // compact arrays:
        uint32_t iw = 0; // write index
        for (uint32_t ir=0; ir<N; ++ir) { // ir = read index
            if (alive[ir]) {
                if (ir>iw) {
                    std::memcpy(&(*this)(0,iw), &(*this)(0,ir), M*sizeof(T));
//                    for (int d=0; d<M; ++d) {
//                        (*this)(d,iw) = operator()(d,ir);
//                    }
                }
                ++iw;
            }
        }
        resize(M,iw);
    }
    void clear() { M=0; N=0; }
    void swap(Matrix<T>& other) {
        uint32_t Mswap = M;
        uint32_t Nswap = N;
        M = other.M;
        N = other.N;
        other.M = Mswap;
        other.N = Nswap;
        T* swapX = X;
        X = other.X;
        other.X = swapX;
        uint32_t X_size_swap = X_size;
        X_size = other.X_size;
        other.X_size = X_size_swap;
    }
};


#endif /* stats_hpp */
