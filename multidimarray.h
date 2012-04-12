#include <iostream>
#include <cstring>
#include <gsl/gsl_vector.h>

class MultiDimArrayLong
{

  private:
    gsl_vector_int* bins;
    long* array;
    long array_length;
    gsl_vector_int* index_multipliers;
    
    long get_array_index(gsl_vector_int* b);

  public:
    MultiDimArrayLong(gsl_vector_int* b);    
    ~MultiDimArrayLong();
    int dim();
    long get(gsl_vector_int* b);
    void set(gsl_vector_int* b, long value);
    void set_all(long value);
    void clear();
    long total();
    void inc(gsl_vector_int* b, long value=1); // +1
    void dec(gsl_vector_int* b, long value=1); // -1
    void print_debug_info();
    // bool test_if_values_are_in_range(long min, long max);
};
