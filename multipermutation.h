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

// Based on the MultiDimArrayLong class, the MultiPermutation class can
// store numbers of type long given a GSL vector that contain multiple
// permutations. For example, when initialized with two permutations of
// length 2 and 4, then for access one provides a vector of length 6
// with entries like (1, 0, 3, 2, 0, 1).

#include <iostream>
#include <cstring>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_sort_vector_double.h>
#include "multidimarray.h"

#define __MULTIPERMUTATION_H

class MultiPermutation {

  private:
    gsl_vector_int* permutation_elements;
    gsl_vector_int* temp_access_vector;
    gsl_vector_int* temp_access_vector_reduced;
    MultiDimArrayLong* mem;
    
    int required_length_of_access_vector() const;
    int required_length_of_reduced_access_vector() const;

    void set_temp_vector_to_upper_bound_of_permutation_values();
    void compute_single_permutation_element(const gsl_vector* input_vector, int in_start, int in_end, gsl_vector_int* resulting_ranks, int res_start=0);
    void set_reduced_temp_vector_to_reduced_upper_bound_of_permutation_values();
    void set_reduced_temp_vector_to_reduced_access_vector(const gsl_vector_int* access);

  public:
    MultiPermutation(const gsl_vector_int* ps);    
    ~MultiPermutation();
    long get(const gsl_vector_int* access, bool assume_validity_of_access=false);
    void set(const gsl_vector_int* access, long value, bool assume_validity_of_access=false);
    // void set_all(long value);
    void clear();
    long total();
    void inc(const gsl_vector_int* access, long value=1, bool assume_validity_of_access=false);
    void dec(const gsl_vector_int* access, long value=1, bool assume_validity_of_access=false);
    void print_debug_info();
    bool test_validity_of_given_access_vector(const gsl_vector_int* access);
    void compute_permutations(const gsl_vector* input_vector, gsl_vector_int* resulting_ranks);
    void write_upper_bound_of_permutation_values_to_vector(gsl_vector_int* output);
};

// General utility functions
void compute_permutation(const gsl_vector* vector, gsl_vector_int* resulting_ranks, int start_index=0);
void compute_permutation(const gsl_vector* vector, gsl_permutation* vector_sorting, gsl_vector_int* resulting_ranks, int start_index=0);
