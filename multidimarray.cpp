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

MultiDimArrayLong::MultiDimArrayLong(const gsl_vector_int* b)
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

long MultiDimArrayLong::get(const gsl_vector_int* b) const
{
  return array[get_array_index(b)];
}

void MultiDimArrayLong::set(const gsl_vector_int* b, long value)
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

long MultiDimArrayLong::total() const
{
  long sum = 0;
  for(long l=0; l<array_length; l++) {
    sum += array[l];
  }
  return sum;
}

void MultiDimArrayLong::inc(const gsl_vector_int* b, long value){
  array[get_array_index(b)] += value;
}
void MultiDimArrayLong::dec(const gsl_vector_int* b, long value){
  array[get_array_index(b)] -= value;
}

void MultiDimArrayLong::print_debug_info() const {
  std::cout <<" ------ DEBUG info: MultiDimArrayLong ------"<<std::endl;
  
  std::cout <<"array_length = "<<array_length<<std::endl;
  std::cout <<"length of bins = "<<bins->size<<std::endl;

  std::cout <<"index_multipliers = ( ";
  for(int i=0; i<index_multipliers->size; i++) {
    std::cout <<gsl_vector_int_get(index_multipliers,i)<<" ";
  }
  std::cout <<")"<<std::endl;
  
  std::cout <<"array = ( ";
  for(long i=0; i<array_length; i++) {
    std::cout <<array[i]<<" ";
  }
  std::cout <<")"<<std::endl;
}

int MultiDimArrayLong::dim() const {
  return bins->size;
}

long MultiDimArrayLong::get_array_index(const gsl_vector_int* b) const
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
      std::cout <<"error: access vector element ("<<gsl_vector_int_get(b,i)<<") at position "<<i<<" exceeds or equals MultiDimArray bounds ("<<gsl_vector_int_get(bins,i)<<"), a.k.a. segmentation fault!"<<std::endl;
      exit(1);
    }
    tempindex += gsl_vector_int_get(b,i)*gsl_vector_int_get(index_multipliers,i);
  }
  // std::cout <<"debug in get_array_index: return value is "<<tempindex<<std::endl;
  return tempindex;
}

long MultiDimArrayLong::get_raw_array_length() {
  return array_length;
}

long MultiDimArrayLong::memory_usage_in_bytes() {
  return array_length * long(sizeof(long));
}

long& MultiDimArrayLong::operator[](const gsl_vector_int* b)
{
  return array[get_array_index(b)];
}
