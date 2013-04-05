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
#include "../multipermutation.h"
#include "../te-datainit.h"
#include "../miniprofiler.h"

using namespace std;


// ----------------------------------------- Tests for MultiDimArrayLong class -----------------------------------------

TEST_CASE("multidimarray/constructors", "Should initialize and construct MultiDimArrayLong object.") {
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

TEST_CASE("multidimarray/input_output", "Should read and write MultiDimArrayLong object.") {
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


TEST_CASE("multidimarray/bulk_methods", "Should bulk edit MultiDimArrayLong object.") {
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

TEST_CASE("multidimarray/increase_and_decrease", "Should edit MultiDimArrayLong entries relatively.") {
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

TEST_CASE("multidimarray/constructor_minimal", "Should work with minimal MultiDimArrayLong object.") {
  gsl_vector_int * access = gsl_vector_int_alloc(1);
  gsl_vector_int_set(access,0,1);
  MultiDimArrayLong* tensor = new MultiDimArrayLong(access);

  CHECK(tensor->dim() == 1);
  CHECK(tensor->total() == 0);
  
  gsl_vector_int_set(access,0,0);
  tensor->inc(access,4);
  CHECK(tensor->total() == 4);
}

TEST_CASE("multidimarray/no_bleeding", "Should set all elements individually and independently.") {
  gsl_vector_int* constructor = gsl_vector_int_alloc(2);
  gsl_vector_int_set(constructor, 0, 2);
  gsl_vector_int_set(constructor, 1, 3);
  MultiDimArrayLong* tensor = new MultiDimArrayLong(constructor);
  tensor->clear();
  
  gsl_vector_int* access = gsl_vector_int_alloc(2);
  // gsl_vector_int_set_zero(access);
  long value = -1;
  // go through all possible indices
  for(int i1=0; i1<2; i1++) {
    gsl_vector_int_set(access, 0, i1);
    for(int j1=0; j1<3; j1++) {
      gsl_vector_int_set(access, 1, j1);
      // test
      CHECK( tensor->get(access) == 0 );
      tensor->set(access, value);
    }
  }
  CHECK( tensor->total() == 2*3*value );
}


// ----------------------------------------- Tests for MultiPermutation class -----------------------------------------

TEST_CASE("multipermutation/constructors", "Should initialize and construct MultiPermutation object.") {
  MultiPermutation* perm = NULL;
  gsl_vector_int * access = NULL;

  access = gsl_vector_int_alloc(1);
  gsl_vector_int_set(access,0,3);
  perm = new MultiPermutation(access);
  perm->clear();

  CHECK(perm->total() == 0);  
}

TEST_CASE("multipermutation/valid_permutations", "Should confirm validity of good MultiPermutation access vectors.") {
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

TEST_CASE("multipermutation/invalid_permutations", "Should reject validity of bad MultiPermutation access vectors.") {
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

// Fails at the moment:
// TEST_CASE("multipermutation/minimal","Should work even with minimal MultiPermutation object.") {
//   MultiPermutation* perm = NULL;
//   gsl_vector_int * access = NULL;
// 
//   access = gsl_vector_int_alloc(1);
//   gsl_vector_int_set(access,0,1);
//   perm = new MultiPermutation(access);
//   perm->clear();
// 
//   CHECK(perm->total() == 0);
//   
//   gsl_vector_int_set(access,0,0);
//   CHECK( perm->test_validity_of_given_access_vector(access) == true );
// 
//   gsl_vector_int_set(access,0,1);
//   CHECK( perm->test_validity_of_given_access_vector(access) == false );
// }

TEST_CASE("multipermutation/set_and_get", "Should set and also return an entry in the MultiPermutation object.") {
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

TEST_CASE("multipermutation/set_and_get_special", "Should set and also return an entry in the MultiPermutation object also in the special case when the first entries of a permutation are the lowest.") {
  MultiPermutation* perm = NULL;
  gsl_vector_int * access = NULL;

  access = gsl_vector_int_alloc(3);
  gsl_vector_int_set(access,0,4);
  gsl_vector_int_set(access,1,4);
  gsl_vector_int_set(access,2,4);
  perm = new MultiPermutation(access);
  
  gsl_vector_int * access_test = gsl_vector_int_alloc(4+4+4);
  gsl_vector_int_set(access_test,0,1);
  gsl_vector_int_set(access_test,1,0);
  gsl_vector_int_set(access_test,2,3);
  gsl_vector_int_set(access_test,3,2);
  gsl_vector_int_set(access_test,4,0);
  gsl_vector_int_set(access_test,5,3);
  gsl_vector_int_set(access_test,6,2);
  gsl_vector_int_set(access_test,7,1);
  gsl_vector_int_set(access_test,8,0);
  gsl_vector_int_set(access_test,9,1);
  gsl_vector_int_set(access_test,10,3);
  gsl_vector_int_set(access_test,11,2);
  long value = 70023;
  perm->set(access_test, value);
  CHECK( perm->get(access_test) == value );
}

TEST_CASE("multipermutation/set_and_get_tolerance", "Should set/get an entry in the MultiPermutation object, even if the access vector is longer than neseccary.") {
  MultiPermutation* perm = NULL;
  gsl_vector_int * access = NULL;

  access = gsl_vector_int_alloc(2);
  gsl_vector_int_set(access,0,2);
  gsl_vector_int_set(access,1,2);
  perm = new MultiPermutation(access);
  
  gsl_vector_int * access_test = gsl_vector_int_alloc(4+3);
  gsl_vector_int_set(access_test,0,1);
  gsl_vector_int_set(access_test,1,0);
  gsl_vector_int_set(access_test,2,0);
  gsl_vector_int_set(access_test,3,1);
  gsl_vector_int_set(access_test,4,2);
  gsl_vector_int_set(access_test,5,0);
  gsl_vector_int_set(access_test,6,-20);
  long value = -34523;
  perm->set(access_test, value);
  CHECK( perm->get(access_test) == value );
}

TEST_CASE("multipermutation/no_bleeding", "Should set all elements individually and independently.") {
  gsl_vector_int* constructor = gsl_vector_int_alloc(2);
  gsl_vector_int_set(constructor, 0, 2);
  gsl_vector_int_set(constructor, 1, 3);
  MultiPermutation* perm = new MultiPermutation(constructor);
  perm->clear();
  
  // cout <<endl<<"DEBUG: --- init ---"<<endl;
  // perm->print_debug_info();
  
  gsl_vector_int* access = gsl_vector_int_alloc(2+3);
  gsl_vector_int_set_zero(access);
  long value = -1;
  // go through all possible permutations
  for(int i1=0; i1<2; i1++) {
    gsl_vector_int_set(access, 0, i1);
    for(int i2=0; i2<2; i2++) {
      gsl_vector_int_set(access, 1, i2);
      for(int j1=0; j1<3; j1++) {
        gsl_vector_int_set(access, 2, j1);
        for(int j2=0; j2<3; j2++) {
          gsl_vector_int_set(access, 3, j2);
          for(int j3=0; j3<3; j3++) {
            gsl_vector_int_set(access, 4, j3);
            // if this is a valid permutation
            if(perm->test_validity_of_given_access_vector(access)) {
              // cout <<endl<<"DEBUG: --- found valid permutation ---"<<endl;
              // cout <<"DEBUG: (i1,i2;j1,j2,j3) = ("<<i1<<","<<i2<<";"<<j1<<","<<j2<<","<<j3<<") => get = "<<perm->get(access)<<endl;
              CHECK( perm->get(access) == 0 );
              perm->set(access, value);
              // perm->print_debug_info();
            }
          }
        }
      }
    }
  }
  CHECK( perm->total() == (2*1)*(3*2*1)*value );
}

TEST_CASE("multipermutation/utility", "Should put the permutation of a given vector into a different one.") {
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

TEST_CASE("multipermutation/utility_custom", "Should work as multipermutation/utility with more options.") {
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

TEST_CASE("multipermutation/compute_permutations", "Should compute multiple permutations.") {
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

TEST_CASE("multipermutation/write_upper_bound_of_permutation_values_to_vector", "Should write limit of perm. values to a vector.") {
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

TEST_CASE("te-datainit/discretize", "Should discretize a continuous signal (without range).") {
  const int len = 5;
  double in[] = {0.1, 1.1, 2.1, 3.1, 4.1};
  rawdata * out = new rawdata[len];
  discretize(&in[0],out,len,len);
  for(int i = 0; i<len; i++) {
    CHECK(int(out[i]) == i);
  }
  
  delete[] out;
}

TEST_CASE("te-datainit/discretize_in_range", "Should discretize a continuous signal (in a given range).") {
  const int len = 5;
  double in[] = {0.1, 1.1, 2.1, 3.1, 4.1};
  rawdata * out = new rawdata[len];
  discretize(&in[0],out,1.1,2.1,len,2);
  for(int i = 0; i<len; i++) {
    CHECK(int(out[i]) == int(in[i]>1.5));
  }
    
  delete[] out;
}

TEST_CASE("te-datainit/discretize_according_to_bin_edges", "Should discretize a continuous signal, given a vector of threshold values.") {
  const int len = 5;
  double in[] = {0.1, 1.1, 2.1, 3.1, 4.1};
  std::vector<double> binEdges; // = {0.0, 1.5, 10000.0}; does not work before C++11
  binEdges.push_back(0.0);
  binEdges.push_back(1.5);
  binEdges.push_back(10000.0);
  
  rawdata * out = new rawdata[len];
  discretize(&in[0], out, len, binEdges);
  for(int i = 0; i<len; i++) {
    CHECK(int(out[i]) == int(in[i]>1.5));
  }
    
  delete[] out;
}

TEST_CASE("te-datainit/compare_discretize_methods", "Should discretize a continuous signal in the same way, given either a vector or just the bin number.") {
  // Goal here: discretize the interval between 1.0 and 2.0 into three bins with two different methods, see if the result matches
  const int len = 14;
  double in[] = {0.1, 1.1, 2.1, 3.1, 4.1, 0.11, 1.11, 2.11, 3.11, 4.11, -0.11, 1.3, 1.333, 1.3334};
  std::vector<double> binEdges; // = {0.0, 1.5, 10000.0}; does not work before C++11
  binEdges.push_back(-1000.0);
  binEdges.push_back(1.333333);
  binEdges.push_back(1.666666);
  binEdges.push_back(1000.0);
  
  rawdata * out_using_bins = new rawdata[len];
  rawdata * out_using_bin_edges = new rawdata[len];
  
  discretize(&in[0], out_using_bins, 1.0, 2.0, len, 3);
  discretize(&in[0], out_using_bin_edges, len, binEdges);
  for(int i = 0; i<len; i++) {
    // cout <<"DEBUG: i = "<<i<<", double value = "<<in[i]<<": out_using_bins = "<<int(out_using_bins[i])<<", out_using_bin_edges = "<<int(out_using_bin_edges[i])<<endl;
    CHECK( out_using_bins[i] == out_using_bin_edges[i] );
  }
  
  delete[] out_using_bins;
  delete[] out_using_bin_edges;
}


TEST_CASE("te-datainit/high_pass", "Should apply a high pass filter to a time series.") {
  const int len = 5;
  double in[] = {0.1, 1.1, 20.1, -30.1, 111.4};
  const double expected[] = {0.0, 1.0, 19.0, -50.2, 141.5};
  const double tolerance = 0.0001;
  
  apply_high_pass_filter_to_time_series(&in[0],len);
  for(int i = 0; i<len; i++) {
    CHECK(abs(in[i]-expected[i])<tolerance);
  }
}


#ifdef SIM_IO_H
TEST_CASE("te-datainit-SKIPPED", "Some of the te-datainit function tests requires the version without the SimKernel bindings!") {
  std::cout <<"Warning: Skipped a number of tests for te-datainit library."<<std::endl;
}
#else

TEST_CASE("te-datainit/generate_from_spikes", "Should generate spike count time series based on spike data.") {
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

TEST_CASE("te-datainit/generate_flouro_from_spikes", "Should generate fluorescence time series based on spike data.") {
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


// ----------------------------------------- Tests for MiniProfiler class -----------------------------------------

TEST_CASE("miniprofiler/constructor", "Should initialize MiniProfiler object.") {
  MiniProfiler prof = MiniProfiler();
  CHECK( prof.number_of_registered_tasks() == 0 );
}

TEST_CASE("miniprofiler/add_task", "Should register new tasks to MiniProfiler object.") {
  MiniProfiler prof = MiniProfiler();
  prof.register_task("some task");
  prof.register_task("another task");
  CHECK( prof.number_of_registered_tasks() == 2 );
}

TEST_CASE("miniprofiler/init_time_zero", "Should register new tasks to MiniProfiler object.") {
  MiniProfiler prof = MiniProfiler();
  prof.register_task("some task");
  CHECK( prof.get_current_time("some task") == 0.0 );
  CHECK( prof.get_current_time("some task that was not even registered") == 0.0 );
}

TEST_CASE("miniprofiler/start_counting", "Should take time of tasks to MiniProfiler object.") {
  MiniProfiler prof = MiniProfiler();
  const std::string label = "some task";
  prof.register_task(label);
  prof.resuming_task(label);
  sleep(0.2);
  prof.stopping_task(label);
  CHECK( prof.get_current_time(label) > 0.0 );
  // cout <<"DEBUG: "<<prof.summary()<<endl;
}

TEST_CASE("miniprofiler/current_time", "Should include time after last stopping of tasks to MiniProfiler object.") {
  MiniProfiler prof = MiniProfiler();
  const std::string label = "some task";
  prof.register_task(label);
  prof.resuming_task(label);
  sleep(0.2);
  prof.stopping_task(label);
  float first_time = prof.get_current_time(label);
  prof.resuming_task(label);
  sleep(0.2);
  CHECK( prof.get_current_time(label) > first_time );
}
