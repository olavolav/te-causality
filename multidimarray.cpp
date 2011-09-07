#include "multidimarray.h"

MultiDimArrayLong::MultiDimArrayLong(gsl_vector_int* b)
// : bins(b)
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
  std::cout <<"debug: allocating array of length "<<supposed_length<<"..."<<std::endl;
  
  // try to allocate memory
  try {
    array = new long[supposed_length];
    memset(array,0,supposed_length*sizeof(double));
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

void MultiDimArrayLong::inc(gsl_vector_int* b, long value){
  array[get_array_index(b)] += value;
}
void MultiDimArrayLong::dec(gsl_vector_int* b, long value){
  array[get_array_index(b)] -= value;
}

void MultiDimArrayLong::dump_debug_info() {
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

// void MultiDimArrayLong::print_slice(int d1=0, int d2=1)
// {
//   if(d1<0 || d1>=dims->size || d2<0 || d2>=dims->size) {
//     std::cout <<"warning in print_slice: invalid dimension settings! skipping function."<<std::endl;
//   }
//   else {
//     for(int i=0; i<4; i++) {
//       gsl_vector_int_set(coordinate,1,i);
//       for(int j=0; j<4; j++) {
//         // inverted to j->0 because we want the first axis as the rows in the output
//         gsl_vector_int_set(coordinate,0,j);
//         cout <<myTensor->get(coordinate)<<"\t";
//       }
//       cout <<endl;
//     }
//   }
// }

int MultiDimArrayLong::dim() {
  return bins->size;
}

long MultiDimArrayLong::get_array_index(gsl_vector_int* b)
{
  long tempindex = 0;
  if(b->size != bins->size) {
    std::cout <<"error: dimensionality of MultiDimArray does not match the one of access vector!"<<std::endl;
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
