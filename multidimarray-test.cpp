// test with:
// 1.) g++ -c multidimarray.cpp -o mda.o
// 2.) g++ multidimarray-test.cpp mda.o -o multidimarray-test -lgsl
// 3.) ./multidimarray-test

#include <iostream>
#include <gsl/gsl_vector.h>
#include "multidimarray.h"

using namespace std;

int main(void)
{
  cout <<"--- testing MultiDimArray object ---"<<endl;
  
  gsl_vector_int * dims = gsl_vector_int_alloc(3);
  gsl_vector_int_set(dims,0,4);
  gsl_vector_int_set(dims,1,3);
  gsl_vector_int_set(dims,2,2);
  MultiDimArrayLong* myTensor = new MultiDimArrayLong(dims);
  // MultiDimArrayLong myTensor(dims);
  
  gsl_vector_int * coordinate = gsl_vector_int_calloc(3);  
  
  cout <<"pure matrix after init:"<<endl;
  for(int i=0; i<3; i++) {
    gsl_vector_int_set(coordinate,1,i);
    for(int j=0; j<4; j++) {
      gsl_vector_int_set(coordinate,0,j);
      gsl_vector_int_set(coordinate,2,0);
      cout <<myTensor->get(coordinate)<<",";
      gsl_vector_int_set(coordinate,2,1);
      cout <<myTensor->get(coordinate)<<"\t";
    }
    cout <<endl;
  }
  myTensor->dump_debug_info();
  
  cout <<"pure matrix after change:"<<endl;
  gsl_vector_int_set_all(coordinate,1);
  // gsl_vector_int_set(coordinate,1,2);
  myTensor->inc(coordinate);
  
  // for(int i=0; i<3; i++) {
  //   gsl_vector_int_set(coordinate,1,i);
  //   for(int j=0; j<4; j++) {
  //     // inverted to j->0 because we want the first axis as the rows in the output
  //     gsl_vector_int_set(coordinate,0,j);
  //     cout <<myTensor->get(coordinate)<<"\t";
  //   }
  //   cout <<endl;
  // }
  // myTensor->dump_debug_info();
  for(int i=0; i<3; i++) {
    gsl_vector_int_set(coordinate,1,i);
    for(int j=0; j<4; j++) {
      gsl_vector_int_set(coordinate,0,j);
      gsl_vector_int_set(coordinate,2,0);
      cout <<myTensor->get(coordinate)<<",";
      gsl_vector_int_set(coordinate,2,1);
      cout <<myTensor->get(coordinate)<<"\t";
    }
    cout <<endl;
  }
  
  cout <<"--- done. ---"<<endl;
  return 0;
}