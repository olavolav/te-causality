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
    
     // copy constructor to duplicate the object (deep copy)
    MultiDimArrayLong(MultiDimArrayLong& copy_from_me);
    long* get_raw_data();
    gsl_vector_int* get_raw_bins_vector();
    long get_raw_array_length();
};
