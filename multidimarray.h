#include <iostream>
#include <gsl/gsl_vector.h>

class MultiDimArrayLong
{
  public:
    // MultiDimArrayLong();    
    MultiDimArrayLong(gsl_vector_int* b);    
    ~MultiDimArrayLong();
    int dim();
    long get(gsl_vector_int* b);
    void set(gsl_vector_int* b, long value);
    void set_all(long value);
    void clear();
    void inc(gsl_vector_int* b, long value=1); // +1
    void dec(gsl_vector_int* b, long value=1); // -1
    void dump_debug_info();
    // void print_slice(int d1=0, int d2=1);
  
  private:
    gsl_vector_int* bins;
    long* array;
    long array_length;
    gsl_vector_int* index_multipliers;
    // long tempindex;
    long get_array_index(gsl_vector_int* b);
    
};
