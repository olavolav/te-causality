// Copyright 2012, Olav Stetter
// 
// This file is part of TE-Causality.
// 
// TE-Causality is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// TE-Causality is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with TE-Causality.  If not, see <http://www.gnu.org/licenses/>.

// This class can store numbers of type long just as long[][]... The
// difference is that the elements are accessed using a GSL vector and
// thus the dimensionality and length of individual arrays can be set
// dynamically. Internally, the vector of indices is mapped to a single
// number, so it uses a one-dimensional array internally.

#include <iostream>
#include <cstring>
#include <gsl/gsl_vector.h>

#define __MULTIDIMARRAYLONG_H

class MultiDimArrayLong {

  private:
    gsl_vector_int* bins;
    long* array;
    long array_length;
    gsl_vector_int* index_multipliers;
    
    long get_array_index(const gsl_vector_int* b) const;

  public:
    MultiDimArrayLong(const gsl_vector_int* b);    
    ~MultiDimArrayLong();
    int dim() const;
    long get(const gsl_vector_int* b) const;
    void set(const gsl_vector_int* b, long value);
    void set_all(long value);
    void clear();
    long total() const;
    void inc(const gsl_vector_int* b, long value=1); // +1
    void dec(const gsl_vector_int* b, long value=1); // -1
    void print_debug_info() const;
    long get_raw_array_length();
    long memory_usage_in_bytes();
    // const long operator[](const gsl_vector_int* b) const;
    long& operator[](const gsl_vector_int* b);
};
