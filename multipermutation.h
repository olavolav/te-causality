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
    MultiDimArrayLong* mem;
    
    long get_array_index(gsl_vector_int* ps);
    int required_length_of_access_vector();
    // int factorial(int n);
    // int encode_single_permutation_to_local_index(gsl_vector_int* perm, int first_index, int last_index);
    // long maximum_value_for_given_index_based_on_access_vector(int dim);
    void set_temp_vector_to_upper_bound_of_permutation_values();
    void compute_single_permutation_element(gsl_vector* input_vector, int in_start, int in_end, gsl_vector_int* resulting_ranks, int res_start);

  public:
    MultiPermutation(gsl_vector_int* ps);    
    ~MultiPermutation();
    long get(gsl_vector_int* access);
    void set(gsl_vector_int* access, long value);
    // void set_all(long value);
    void clear();
    long total();
    void inc(gsl_vector_int* access, long value=1); // +1
    void dec(gsl_vector_int* access, long value=1); // -1
    void print_debug_info();
    bool test_validity_of_given_access_vector(gsl_vector_int* access);
    void compute_permutations(gsl_vector* input_vector, gsl_vector_int* resulting_ranks);
    void write_upper_bound_of_permutation_values_to_vector(gsl_vector_int* output);
};

// General utility functions
void compute_permutation(gsl_vector* vector, gsl_vector_int* resulting_ranks, int start_index=0);
void compute_permutation(gsl_vector* vector, gsl_permutation* vector_sorting, gsl_vector_int* resulting_ranks, int start_index=0);
