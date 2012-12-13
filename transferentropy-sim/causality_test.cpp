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

#define CATCH_CONFIG_MAIN  // This tell CATCH to provide a main() - only do this in one cpp file
#include "tests/catch.hpp"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
// #include "../multidimarray.h"
#include "../multipermutation.h"
#include "../te-datainit.h"

using namespace std;


// ----------------------------------------- Tests for MultiDimArrayLong class -----------------------------------------

TEST_CASE("multidimarray/constructors","Should initialize and construct MultiDimArrayLong object.") {
  MultiDimArrayLong* tensor = NULL;
  gsl_vector_int * access = NULL;
  
  access = gsl_vector_int_alloc(3);
  gsl_vector_int_set(access,0,3);
  gsl_vector_int_set(access,1,4);
  gsl_vector_int_set(access,2,5);
  tensor = new MultiDimArrayLong(access);

  CHECK(tensor->dim() == 3);
  CHECK(tensor->total() == 0);  
}

TEST_CASE("multidimarray/input_output","Should read and write MultiDimArrayLong object.") {
  MultiDimArrayLong* tensor = NULL;
  gsl_vector_int * access = NULL;
  
  access = gsl_vector_int_alloc(3);
  gsl_vector_int_set(access,0,3);
  gsl_vector_int_set(access,1,4);
  gsl_vector_int_set(access,2,5);
  tensor = new MultiDimArrayLong(access);

  for(int x=0;x<3;x++) {
    gsl_vector_int_set(access,0,x);
    for(int y=0;y<4;y++) {
      gsl_vector_int_set(access,1,y);
      for(int z=0;z<5;z++) {
        gsl_vector_int_set(access,2,z);
        tensor->set(access,(x+1)*(7*y+1)*(11*z+1)%13);
      }
    }
  }
  for(int x=0;x<3;x++) {
    gsl_vector_int_set(access,0,x);
    for(int y=0;y<4;y++) {
      gsl_vector_int_set(access,1,y);
      for(int z=0;z<5;z++) {
        gsl_vector_int_set(access,2,z);
        CHECK(tensor->get(access) == (x+1)*(7*y+1)*(11*z+1)%13);
      }
    }
  }
}


TEST_CASE("multidimarray/bulk_methods","Should bulk edit MultiDimArrayLong object.") {
  MultiDimArrayLong* tensor = NULL;
  gsl_vector_int * access = NULL;
  
  access = gsl_vector_int_alloc(3);
  gsl_vector_int_set(access,0,3);
  gsl_vector_int_set(access,1,4);
  gsl_vector_int_set(access,2,5);
  tensor = new MultiDimArrayLong(access);

  tensor->clear();
  CHECK(tensor->total() == 0);
}

TEST_CASE("multidimarray/increase_and_decrease","Should edit MultiDimArrayLong entries relatively.") {
  MultiDimArrayLong* tensor = NULL;
  gsl_vector_int * access = NULL;
  
  access = gsl_vector_int_alloc(2);
  gsl_vector_int_set(access,0,2);
  gsl_vector_int_set(access,1,4);
  tensor = new MultiDimArrayLong(access);
  tensor->clear();

  gsl_vector_int_set(access,0,0);
  gsl_vector_int_set(access,1,2);
  tensor->inc(access);
  CHECK( tensor->get(access) == 1 );
  CHECK( tensor->total() == 1 );

  gsl_vector_int_set(access,0,1);
  gsl_vector_int_set(access,1,0);
  tensor->dec(access, 4);
  CHECK( tensor->total() == -3 );
}


// ----------------------------------------- Tests for MultiPermutation class -----------------------------------------

TEST_CASE("multipermutation/constructors","Should initialize and construct MultiPermutation object.") {
  MultiPermutation* perm = NULL;
  gsl_vector_int * access = NULL;

  access = gsl_vector_int_alloc(1);
  gsl_vector_int_set(access,0,3);
  perm = new MultiPermutation(access);
  perm->clear();

  CHECK(perm->total() == 0);  
}

TEST_CASE("multipermutation/valid_permutations","Should confirm validity of good MultiPermutation access vectors.") {
  MultiPermutation* perm = NULL;
  gsl_vector_int * access = NULL;

  access = gsl_vector_int_alloc(3);
  gsl_vector_int_set(access,0,2);
  gsl_vector_int_set(access,1,1);
  gsl_vector_int_set(access,2,5);
  perm = new MultiPermutation(access);
  
  gsl_vector_int * access_test = gsl_vector_int_alloc(2+1+5);
  gsl_vector_int_set(access_test,0,1);
  gsl_vector_int_set(access_test,1,0);
  gsl_vector_int_set(access_test,2,0);
  gsl_vector_int_set(access_test,3,2);
  gsl_vector_int_set(access_test,4,0);
  gsl_vector_int_set(access_test,5,1);
  gsl_vector_int_set(access_test,6,4);
  gsl_vector_int_set(access_test,7,3);
  CHECK( perm->test_validity_of_given_access_vector(access_test) == true );

  gsl_vector_int_set(access_test,0,0);
  gsl_vector_int_set(access_test,1,1);
  // gsl_vector_int_set(access_test,2,0);
  gsl_vector_int_set(access_test,3,4);
  gsl_vector_int_set(access_test,4,3);
  gsl_vector_int_set(access_test,5,2);
  gsl_vector_int_set(access_test,6,1);
  gsl_vector_int_set(access_test,7,0);
  CHECK( perm->test_validity_of_given_access_vector(access_test) == true );
}

TEST_CASE("multipermutation/invalid_permutations","Should reject validity of bad MultiPermutation access vectors.") {
  MultiPermutation* perm = NULL;
  gsl_vector_int * access = NULL;

  access = gsl_vector_int_alloc(3);
  gsl_vector_int_set(access,0,2);
  gsl_vector_int_set(access,1,1);
  gsl_vector_int_set(access,2,5);
  perm = new MultiPermutation(access);
  
  gsl_vector_int * access_test = gsl_vector_int_alloc(2+1+5);
  gsl_vector_int_set(access_test,0,1);
  gsl_vector_int_set(access_test,1,0);
  gsl_vector_int_set(access_test,2,1);
  gsl_vector_int_set(access_test,3,2);
  gsl_vector_int_set(access_test,4,0);
  gsl_vector_int_set(access_test,5,1);
  gsl_vector_int_set(access_test,6,4);
  gsl_vector_int_set(access_test,7,3);
  CHECK( perm->test_validity_of_given_access_vector(access_test) == false );

  gsl_vector_int_set(access_test,0,0);
  gsl_vector_int_set(access_test,1,1);
  gsl_vector_int_set(access_test,2,0);
  gsl_vector_int_set(access_test,3,4);
  gsl_vector_int_set(access_test,4,3);
  gsl_vector_int_set(access_test,5,2);
  gsl_vector_int_set(access_test,6,2);
  gsl_vector_int_set(access_test,7,0);
  CHECK( perm->test_validity_of_given_access_vector(access_test) == false );
}

TEST_CASE("multipermutation/minimal","Should work even with minimal MultiPermutation object.") {
  MultiPermutation* perm = NULL;
  gsl_vector_int * access = NULL;

  access = gsl_vector_int_alloc(1);
  gsl_vector_int_set(access,0,1);
  perm = new MultiPermutation(access);
  perm->clear();

  CHECK(perm->total() == 0);
  
  gsl_vector_int_set(access,0,0);
  CHECK( perm->test_validity_of_given_access_vector(access) == true );

  gsl_vector_int_set(access,0,1);
  CHECK( perm->test_validity_of_given_access_vector(access) == false );
}

TEST_CASE("multipermutation/set_and_get","Should set and also return an entry in the MultiPermutation object.") {
  MultiPermutation* perm = NULL;
  gsl_vector_int * access = NULL;

  access = gsl_vector_int_alloc(2);
  gsl_vector_int_set(access,0,2);
  gsl_vector_int_set(access,1,2);
  perm = new MultiPermutation(access);
  
  gsl_vector_int * access_test = gsl_vector_int_alloc(2+2);
  gsl_vector_int_set(access_test,0,1);
  gsl_vector_int_set(access_test,1,0);
  gsl_vector_int_set(access_test,2,0);
  gsl_vector_int_set(access_test,3,1);
  long value = 123134523;
  perm->set(access_test, value);
  CHECK( perm->get(access_test) == value );
}

TEST_CASE("multipermutation/utility","Should put the permutation of a given vector into a different one.") {
  gsl_vector* time_series = gsl_vector_alloc(4);
  gsl_vector_set(time_series, 0, 1.01);
  gsl_vector_set(time_series, 1, 3.14);
  gsl_vector_set(time_series, 2, 2.75);
  gsl_vector_set(time_series, 3, 1.02);

  gsl_vector_int* ranks = gsl_vector_int_alloc(4);
  gsl_vector_int_set_all(ranks, -1);
  compute_permutation(time_series,ranks);
  
  CHECK( gsl_vector_int_get(ranks, 0) == 0 );
  CHECK( gsl_vector_int_get(ranks, 1) == 3 );
  CHECK( gsl_vector_int_get(ranks, 2) == 2 );
  CHECK( gsl_vector_int_get(ranks, 3) == 1 );
}

TEST_CASE("multipermutation/utility_custom","Should work as multipermutation/utility with more options.") {
  gsl_vector* time_series = gsl_vector_alloc(3);
  gsl_vector_set(time_series, 0, 10.01);
  gsl_vector_set(time_series, 1, -3.14);
  gsl_vector_set(time_series, 2, 2.75);

  gsl_vector_int* ranks = gsl_vector_int_alloc(6);
  gsl_vector_int_set_all(ranks, -1);
  compute_permutation(time_series,ranks,1);
  
  CHECK( gsl_vector_int_get(ranks, 0) == -1 );
  CHECK( gsl_vector_int_get(ranks, 1) == 2 );
  CHECK( gsl_vector_int_get(ranks, 2) == 0 );
  CHECK( gsl_vector_int_get(ranks, 3) == 1 );
  CHECK( gsl_vector_int_get(ranks, 4) == -1 );
  CHECK( gsl_vector_int_get(ranks, 5) == -1 );
}

TEST_CASE("multipermutation/compute_permutations","Should compute multiple permutations.") {
  MultiPermutation* perm = NULL;
  gsl_vector_int * access = NULL;

  access = gsl_vector_int_alloc(2);
  gsl_vector_int_set(access,0,2);
  gsl_vector_int_set(access,1,4);
  perm = new MultiPermutation(access);
  const int length = 2+4;

  gsl_vector * time_series = gsl_vector_alloc(length);
  gsl_vector_set(time_series,0,1.23);
  gsl_vector_set(time_series,1,-1.23);
  gsl_vector_set(time_series,2,1.23);
  gsl_vector_set(time_series,3,10.23);
  gsl_vector_set(time_series,4,100.23);
  gsl_vector_set(time_series,5,-1.23);
  
  gsl_vector_int * permutations = gsl_vector_int_alloc(length);
  perm->compute_permutations(time_series, permutations);
  
  // First, check that the result is a valid access vector
  CHECK( perm->test_validity_of_given_access_vector(permutations) );
  
  // Then, make sure the entries are correctly ordered
  gsl_vector_int * access_test = gsl_vector_int_alloc(length);
  CHECK( gsl_vector_int_get(permutations, 0) == 1 );
  CHECK( gsl_vector_int_get(permutations, 1) == 0 );
  CHECK( gsl_vector_int_get(permutations, 2) == 1 );
  CHECK( gsl_vector_int_get(permutations, 3) == 2 );
  CHECK( gsl_vector_int_get(permutations, 4) == 3 );
  CHECK( gsl_vector_int_get(permutations, 5) == 0 );
}

TEST_CASE("multipermutation/write_upper_bound_of_permutation_values_to_vector","Should write limit of perm. values to a vector.") {
  MultiPermutation* perm = NULL;
  gsl_vector_int * generator = NULL;

  generator = gsl_vector_int_alloc(2);
  gsl_vector_int_set(generator,0,3);
  gsl_vector_int_set(generator,1,2);
  perm = new MultiPermutation(generator);
  
  gsl_vector_int * access_test = gsl_vector_int_alloc(3+2);
  perm->write_upper_bound_of_permutation_values_to_vector(access_test);
  CHECK( gsl_vector_int_get(access_test,0) == 3 );
  CHECK( gsl_vector_int_get(access_test,1) == 3 );
  CHECK( gsl_vector_int_get(access_test,2) == 3 );
  CHECK( gsl_vector_int_get(access_test,3) == 2 );
  CHECK( gsl_vector_int_get(access_test,4) == 2 );
}


// ----------------------------------------- Tests for te-datainit library -----------------------------------------

#ifdef SIM_IO_H
TEST_CASE("te-datainit-SKIPPED","For the moment, this tests to choose the te-datainit version without the SimKernel bindings!") {
  std::cout <<"Warning: Skipped tests for te-datainit library."<<std::endl;
  }
#else
TEST_CASE("te-datainit/discretize","Should discretize a continuous signal (without range).") {
  const int len = 5;
  double in[] = {0.1, 1.1, 2.1, 3.1, 4.1};
  rawdata * out = new rawdata[len];
  discretize(&in[0],out,len,len);
  for(int i = 0; i<len; i++) {
    CHECK(int(out[i]) == i);
  }
  
  delete[] out;
}

TEST_CASE("te-datainit/discretize_in_range","Should discretize a continuous signal (in a given range).") {
  const int len = 5;
  double in[] = {0.1, 1.1, 2.1, 3.1, 4.1};
  rawdata * out = new rawdata[len];
  discretize(&in[0],out,1.1,2.1,len,2);
  for(int i = 0; i<len; i++) {
    CHECK(int(out[i]) == int(in[i]>1.5));
  }
    
  delete[] out;
}

TEST_CASE("te-datainit/high_pass","Should apply a high pass filter to a time series.") {
  const int len = 5;
  double in[] = {0.1, 1.1, 20.1, -30.1, 111.4};
  const double expected[] = {0.0, 1.0, 19.0, -50.2, 141.5};
  const double tolerance = 0.0001;
  
  apply_high_pass_filter_to_time_series(&in[0],len);
  for(int i = 0; i<len; i++) {
    CHECK(abs(in[i]-expected[i])<tolerance);
  }
}

TEST_CASE("te-datainit/generate_from_spikes","Should generate spike count time series based on spike data.") {
  string inputfile_times = "tests/spiketimes_test.txt";
  string inputfile_indices = "tests/spikeindices_test.txt";
  int size = 2;
  int tauImg = 5.; // ms
  long samples = 4;
  string model = "SpikeCount"; //"Leogang"
  double noise = -1;
  double saturation = -1.0;
  double cutoff = -1.0;
  double deltaCa = 50.0;
  double tauCa = 100000.0; // a bit slow
  const double tolerance = 0.0001;
  
  // initialize random number generator
  gsl_rng_env_setup();
  gsl_rng* GSLrandom = gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(GSLrandom, 1234);

  double ** time_series = NULL;
  time_series = generate_time_series_from_spike_data(inputfile_times, inputfile_indices, size, tauImg, samples, model, noise, saturation, cutoff, deltaCa, tauCa, GSLrandom);
  // display_subset(time_series[0],samples);
  // display_subset(time_series[1],samples);
  const double expected[] = {0,1,0,0};
  for(int i = 0; i<samples; i++) {
    CHECK(abs(time_series[0][i]-expected[i])<tolerance);
    CHECK(abs(time_series[1][samples-i-1]-expected[i])<tolerance);
    // ... because the second time series should be the first, just mirrored
  }
}

TEST_CASE("te-datainit/generate_flouro_from_spikes","Should generate fluorescence time series based on spike data.") {
  string inputfile_times = "tests/spiketimes_test.txt";
  string inputfile_indices = "tests/spikeindices_test.txt";
  int size = 2;
  int tauImg = 5.; // ms
  long samples = 4;
  string model = "Leogang";
  double noise = 0.01;
  double saturation = -1.0;
  double cutoff = -1.0;
  double deltaCa = 50.0;
  double tauCa = 100000.0; // a bit slow
  const double tolerance = 1.0;
  
  // initialize random number generator
  gsl_rng_env_setup();
  gsl_rng* GSLrandom = gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(GSLrandom, 1234);

  double ** time_series = NULL;
  time_series = generate_time_series_from_spike_data(inputfile_times, inputfile_indices, size, tauImg, samples, model, noise, saturation, cutoff, deltaCa, tauCa, GSLrandom);
  // display_subset(time_series[0],samples);
  // display_subset(time_series[1],samples);
  const double expected[] = {0,50.0,50.0,50.0};
  for(int i = 0; i<samples; i++) {
    CHECK(abs(time_series[0][i]-expected[i])<tolerance);
  }
}

#endif