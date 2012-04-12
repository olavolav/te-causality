#define CATCH_CONFIG_MAIN  // This tell CATCH to provide a main() - only do this in one cpp file
#include "tests/catch.hpp"

#include <gsl/gsl_vector.h>
#include "../multidimarray.h"
#include "../te-datainit.h"

using namespace std;



// ----------------------------------------- Tests for MultiDimArrayLong class -----------------------------------------

TEST_CASE("multidimarray/constructors","Should initiaize and construct MultiDimArrayLong object.") {
  MultiDimArrayLong* tensor = NULL;
  gsl_vector_int * access = NULL;
  
  access = gsl_vector_int_alloc(3);
  gsl_vector_int_set(access,0,3);
  gsl_vector_int_set(access,1,4);
  gsl_vector_int_set(access,2,5);
  tensor = new MultiDimArrayLong(access);

  REQUIRE(tensor->dim() == 3);
  REQUIRE(tensor->total() == 0);  
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
        REQUIRE(tensor->get(access) == (x+1)*(7*y+1)*(11*z+1)%13);
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
  REQUIRE(tensor->total() == 0);
}



// ----------------------------------------- Tests for te-datainit library -----------------------------------------

#ifdef SIM_IO_H
  TEST_CASE("te-datainit-SKIPPED","For the moment, this requires to choose the te-datainit version without the SimKernel bindings!") {
    std::cout <<"Warning: Skipped tests for te-datainit library."<<std::endl;
  }
#else
TEST_CASE("te-datainit/discretize","Should discretize a continuous signal (without range).") {
  const int len = 5;
  double in[] = {0.1, 1.1, 2.1, 3.1, 4.1};
  rawdata * out = new rawdata[len];
  discretize(&in[0],out,len,len);
  for(int i = 0; i<len; i++) {
    REQUIRE(int(out[i]) == i);
  }
    
  delete[] out;
}

TEST_CASE("te-datainit/discretize_in_range","Should discretize a continuous signal (in a given range).") {
  const int len = 5;
  double in[] = {0.1, 1.1, 2.1, 3.1, 4.1};
  rawdata * out = new rawdata[len];
  discretize(&in[0],out,1.1,2.1,len,2);
  for(int i = 0; i<len; i++) {
    REQUIRE(int(out[i]) == int(in[i]>1.5));
  }
    
  delete[] out;
}

#endif