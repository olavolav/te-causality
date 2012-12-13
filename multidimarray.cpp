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

#include "multidimarray.h"

MultiDimArrayLong::MultiDimArrayLong(gsl_vector_int* b)
{
  if(b->size < 1) {
    std::cout <<"error: init vector for MultiDimArrayLong has zero length!"<<std::endl;
    exit(1);
  }
  bins = gsl_vector_int_alloc(b->size);
  gsl_vector_int_memcpy(bins,b);
  
  // determine length of 1-D array to be allocated
  long supposed_length = 1;
  for(int i=0; i<bins->size; i++) {
    if(gsl_vector_int_get(bins,i) < 1) {
      std::cout <<"error: invalid init vector for MultiDimArrayLong!"<<std::endl;
      exit(1);
    }
    supposed_length *= (long)gsl_vector_int_get(bins,i);
  }
  // std::cout <<"debug: allocating array of length "<<supposed_length<<"..."<<std::endl;
  
  // try to allocate memory
  array = NULL;
  try {
    array = new long[supposed_length];
    memset(array,0,supposed_length*sizeof(long));
  }
  catch(...) {
    std::cout <<"error: failed to allocate enough memory for MultiDimArrayLong!"<<std::endl;
    exit(1);
  };
  array_length = supposed_length;
  
  // set up multipliers
  index_multipliers = gsl_vector_int_alloc(bins->size);
  gsl_vector_int_set_all(index_multipliers,1);
  for (int i=1; i<index_multipliers->size; i++) {
    gsl_vector_int_set(index_multipliers,i,gsl_vector_int_get(bins,i-1)*gsl_vector_int_get(index_multipliers,i-1)); 
  }
  
}

MultiDimArrayLong::~MultiDimArrayLong()
{
  // destruct
  if (array != NULL) delete[] array;
  gsl_vector_int_free(bins);
  gsl_vector_int_free(index_multipliers);
}

long MultiDimArrayLong::get(gsl_vector_int* b)
{
  return array[get_array_index(b)];
}

void MultiDimArrayLong::set(gsl_vector_int* b, long value)
{
  array[get_array_index(b)] = value;
}

void MultiDimArrayLong::set_all(long value)
{
  for(long l=0; l<array_length; l++) {
    array[l] = value;
  }
}

void MultiDimArrayLong::clear()
{
  // MultiDimArrayLong::set_all(0);
  memset(array,0,array_length*sizeof(long));
}

long MultiDimArrayLong::total()
{
  long sum = 0;
  for(long l=0; l<array_length; l++) {
    sum += array[l];
  }
  return sum;
}

void MultiDimArrayLong::inc(gsl_vector_int* b, long value){
  array[get_array_index(b)] += value;
}
void MultiDimArrayLong::dec(gsl_vector_int* b, long value){
  array[get_array_index(b)] -= value;
}

void MultiDimArrayLong::print_debug_info() {
  std::cout <<"debug: array_length = "<<array_length<<std::endl;
  std::cout <<"debug: length of bins = "<<bins->size<<std::endl;

  std::cout <<"debug: index_multipliers = ( ";
  for(int i=0; i<index_multipliers->size; i++) {
    std::cout <<gsl_vector_int_get(index_multipliers,i)<<" ";
  }
  std::cout <<")"<<std::endl;
  
  std::cout <<"debug: array = ( ";
  for(long i=0; i<array_length; i++) {
    std::cout <<array[i]<<" ";
  }
  std::cout <<")"<<std::endl;
}

int MultiDimArrayLong::dim() {
  return bins->size;
}

long MultiDimArrayLong::get_array_index(gsl_vector_int* b)
{
  long tempindex = 0;
  // if(b->size != bins->size) {
  //   std::cout <<"error: dimensionality of MultiDimArray does not match the one of access vector!"<<std::endl;
  //   exit(1);
  // }
  if(b->size < bins->size) {
    // We will quietly ignore trailing entries for now, because then we don't need to create new access
    // vectors constantly when working with MultiDimArrays of different dimensionality.
    std::cout <<"error: dimensionality of MultiDimArray is less than the one of access vector!"<<std::endl;
    exit(1);
  }
  
  for(int i=0; i<bins->size; i++) {
    if(gsl_vector_int_get(b,i)<0 || gsl_vector_int_get(b,i)>=gsl_vector_int_get(bins,i)) {
      std::cout <<"error: access vector exceeds MultiDimArray bounds (segmentation fault)!"<<std::endl;
      exit(1);
    }
    tempindex += gsl_vector_int_get(b,i)*gsl_vector_int_get(index_multipliers,i);
  }
  // std::cout <<"debug in get_array_index: return value is "<<tempindex<<std::endl;
  return tempindex;
}


// copy constructor
MultiDimArrayLong::MultiDimArrayLong(MultiDimArrayLong& copy_from_me) {
  MultiDimArrayLong(copy_from_me.get_raw_bins_vector());
  memcpy(array,copy_from_me.get_raw_data(),copy_from_me.get_raw_array_length()*sizeof(long));
}

long* MultiDimArrayLong::get_raw_data() {
  return array;
}

gsl_vector_int* MultiDimArrayLong::get_raw_bins_vector() {
  return bins;
}

long MultiDimArrayLong::get_raw_array_length() {
  return array_length;
}
