// created Do 28 Apr 2011 18:13:44 CEST by olav

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <ctime>
#include <cstring>
#include <sstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_histogram.h>

#define SUBSET_LENGTH 100

typedef unsigned char rawdata;


double** load_time_series_from_binary_file(std::string inputfile_name, unsigned int size, long samples, double input_scaling, bool OverrideRescalingQ, double std_noise, double fluorescence_saturation, double cutoff);

rawdata* generate_discretized_global_time_series(double** time_series, unsigned int size, long samples, unsigned int globalbins, double GlobalConditioningLevel, unsigned long* AvailableSamples, long StartSampleIndex, long EndSampleIndex);

rawdata** generate_discretized_version_of_time_series(double** in, unsigned int size, long nr_samples, unsigned int nr_bins);
// void discretize(double** in, rawdata** out, unsigned int size, long nr_samples, unsigned int nr_bins);
void discretize(double* in, rawdata* out, long nr_samples, unsigned int nr_bins);
void discretize(double* in, rawdata* out, double min, double max, long nr_samples, unsigned int nr_bins);
rawdata discretize(double in, double min, double max, unsigned int nr_bins);
// void discretize2accordingtoStd(double* in, rawdata* out);

void apply_high_pass_filter_to_time_series(double** time_series, unsigned int size, long nr_samples);
void apply_high_pass_filter_to_time_series(double* time_series, long nr_samples);

double** generate_time_series_from_spike_data(std::string inputfile_spiketimes, std::string inputfile_spikeindices, unsigned int size, unsigned int tauImg, long samples, std::string fluorescence_model, double std_noise=-1., double fluorescence_saturation=-1., double cutoff=-1., double DeltaCalciumOnAP=50., double tauCa=1000.);

unsigned long count(int* array, unsigned long starti, unsigned long endi, int what);
bool has_index(int* array, unsigned long starti, unsigned long endi, int what);

double smallest(double* array, long length);
double largest(double* array, long length);
rawdata smallest(rawdata* array, long length);
rawdata largest(rawdata* array, long length);
double smallest(double** array, unsigned int size, long length);
double largest(double** array, unsigned int size, long length);

void free_time_series_memory(double** xresult, unsigned int size);
void free_time_series_memory(double* xresult);
void free_time_series_memory(rawdata** xresult, unsigned int size);
void free_time_series_memory(rawdata* xresult);

void display_subset(double* data);
void display_subset(rawdata* data);

int Magic_GuessBinNumber(double** data, unsigned int size, long samples);
double Magic_GuessConditioningLevel(double** data, unsigned int size, long samples);

// void Test_SetMinimaToZero(double** data, unsigned int size, long samples);
// void Test_PlotLogHistogram(long* histo, int length);
// void Test_PlotHistogram(long* histo, int length);

void PlotHistogramInASCII(double* data, int samples, double range_min, double range_max, std::string xlabel, std::string ylabel);
void PlotLogHistogramInASCII(double* data, int samples, double range_min, double range_max, std::string xlabel, std::string ylabel);
void PlotLogLogHistogramInASCII(double* data, int samples, double range_min, double range_max, std::string xlabel, std::string ylabel);
void PlotHistogramInASCII(bool xlogscaleQ, bool ylogscaleQ, double* data, int samples, double range_min, double range_max, std::string xlabel="unlabeled", std::string ylabel="unlabeled");
