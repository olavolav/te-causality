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
#include <gsl/gsl_fit.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_gamma.h>


#define ENABLE_YAML_IMPORT_AT_COMPILE_TIME

#ifdef ENABLE_YAML_IMPORT_AT_COMPILE_TIME
#include "yaml-cpp/yaml.h"
#endif

// set output stream depending on wether SimKernel's sim.h is included
// (see also te-datainit.cpp)
// Unit tests will only work without SimKernel!
#include <sim_main.h>

#undef IOSTREAMD
#ifdef SIM_IO_H
  #define IOSTREAMD Sim& output
#else
  #define IOSTREAMD std::ostream* output=&(std::cout)
#endif

// constants for the differential entropy
#define EULERGAMMA 0.57721566490153
#define PI 3.1415926535898

// #define SUBSET_LENGTH 30

typedef unsigned char rawdata;


double** load_time_series_from_file(std::string inputfile_name, unsigned int size, long samples, double input_scaling, bool OverrideRescalingQ, double std_noise, double fluorescence_saturation, double cutoff, gsl_rng* GSLrandom, IOSTREAMD);

rawdata* generate_discretized_global_time_series(double** time_series, unsigned int size, long samples, unsigned int globalbins, double GlobalConditioningLevel, unsigned long* AvailableSamples, long StartSampleIndex, long EndSampleIndex, bool EqualSampleNumberQ, long MaxSampleNumberPerBin, IOSTREAMD);

rawdata** generate_discretized_version_of_time_series(double** in, unsigned int size, long nr_samples, unsigned int nr_bins);
// void discretize(double** in, rawdata** out, unsigned int size, long nr_samples, unsigned int nr_bins);
void discretize(double* in, rawdata* out, long nr_samples, unsigned int nr_bins);
void discretize(double* in, rawdata* out, double min, double max, long nr_samples, unsigned int nr_bins);
rawdata discretize(double in, double min, double max, unsigned int nr_bins);
// void discretize2accordingtoStd(double* in, rawdata* out);

// Orlandi: Adding option for predefined binning limits
rawdata** generate_discretized_version_of_time_series(double** in, unsigned int size, long nr_samples, std::vector<double> binEdges);
void discretize(double* in, rawdata* out, long nr_samples, std::vector<double> binEdges);
rawdata discretize(double in, std::vector<double> binEdges);

void apply_high_pass_filter_to_time_series(double** time_series, unsigned int size, long nr_samples);
void apply_high_pass_filter_to_time_series(double* time_series, long nr_samples);

double** generate_time_series_from_spike_data(std::string inputfile_spiketimes, std::string inputfile_spikeindices, unsigned int size, unsigned int tauImg, long samples, std::string fluorescence_model, double std_noise, double fluorescence_saturation, double cutoff, double DeltaCalciumOnAP, double tauCa, gsl_rng* GSLrandom, IOSTREAMD);

unsigned long count(int* array, unsigned long starti, unsigned long endi, int what);
bool has_index(int* array, unsigned long starti, unsigned long endi, int what);

double smallest(double* array, const long length);
double largest(double* array, const long length);
rawdata smallest(rawdata* array, const long length);
rawdata largest(rawdata* array, const long length);
double smallest(double** array, const unsigned int size, const long length);
double largest(double** array, const unsigned int size, const long length);

double total(double* array, const long first, const long last);
double total(double* array, const long length);
double mean(double* array, const long first, const long last);
double mean(double* array, const long length);
double variance(double* array, const long first, const long last);
double variance(double* array, const long length);
double standard_deviation(double* array, const long first, const long last);
double standard_deviation(double* array, const long length);

double* generate_mean_time_series(double** data, unsigned int size, long samples);

void free_time_series_memory(double** xresult, unsigned int size);
void free_time_series_memory(double* xresult);
void free_time_series_memory(rawdata** xresult, unsigned int size);
void free_time_series_memory(rawdata* xresult);

void display_subset(double* data, int length, IOSTREAMD);
void display_subset(rawdata* data, int length, IOSTREAMD);

int Magic_GuessBinNumber(double** data, const unsigned int size, const long samples);
double Magic_GuessConditioningLevel(double** data, unsigned int size, const long samples, IOSTREAMD);

// void Test_SetMinimaToZero(double** data, unsigned int size, long samples);
// void Test_PlotLogHistogram(long* histo, int length);
// void Test_PlotHistogram(long* histo, int length);

void PlotHistogramInASCII(double* data, int samples, double range_min, double range_max, std::string xlabel, std::string ylabel, IOSTREAMD);
void PlotLogHistogramInASCII(double* data, int samples, double range_min, double range_max, std::string xlabel, std::string ylabel, IOSTREAMD);
void PlotLogLogHistogramInASCII(double* data, int samples, double range_min, double range_max, std::string xlabel, std::string ylabel, IOSTREAMD);
void PlotHistogramInASCII(bool xlogscaleQ, bool ylogscaleQ, double* data, int samples, double range_min, double range_max, std::string xlabel, std::string ylabel, IOSTREAMD);

// void PlotInASCII(double* data, int x_length, double range_min, double range_max, IOSTREAMD);
// void PlotInASCII(double** data, int x_length, double range_min, double range_max, IOSTREAMD);

double Util_FindPeakInHistogram(const double* data, const long samples, const double range_min, const double range_max, const int histo_bins);
void Util_CreateFakeLogLogHistogram(double** x, double** y, double** w, const double* data, const long samples, const double range_min, const double range_max, const int histo_bins=40);
void Util_FreeFakeHistogram(double* x, double* y, double* w);

void Util_CoordinatedForMathematica(double*x, double*y, int length, IOSTREAMD);

double** generate_conditioned_time_series_by_glueing(double** data, const int size, double* xmean, const long StartSampleIndex, const long EndSampleIndex, const double condlevel, unsigned long* available_samples, IOSTREAMD);

#ifdef ENABLE_YAML_IMPORT_AT_COMPILE_TIME
double** read_positions_from_YAML(std::string YAMLfilename, unsigned int size, IOSTREAMD);
#endif
void free_position_memory(double** pos, unsigned int size);

double norm(double* pointA, double* pointB);
double** clone_time_series(double** data, unsigned int size, long samples);
void apply_light_scattering_to_time_series(double** data, unsigned int size, long samples, std::string YAMLfilename, double sigma_scatter, double amplitude_scatter, IOSTREAMD);

double SphericalUnitSurface(int r);
double gsl_norm(const gsl_vector* vecA, const gsl_vector* vecB, int dim);
double gsl_quicknorm(const gsl_vector* vecA, const gsl_vector* vecB, int dim, double bound=10^16);
long double DifferentialEntropy(gsl_vector** data, const int dim, const long samples);

long* generate_random_permutation(long samples, gsl_rng* GSLrandom);
long* generate_random_permutation(long samples, rawdata globalbins, unsigned long* AvailableSamples, long StartSampleIndex, long EndSampleIndex, rawdata* xglobal, gsl_rng* GSLrandom);
void random_permutation(long** data, const long samples, gsl_rng* GSLrandom);
void geometric_permutation(long** data, const long samples, const long AutoCorrLength, gsl_rng* GSLrandom);
long* generate_random_geometric_permutation(long samples, rawdata globalbins, rawdata* xglobal, long AutoCorrLength, gsl_rng* GSLrandom);

double AutoCorrelation(double* data, const long samples, const long lag=0, bool Abs=false);
double AutoCorrelationTimeScale(double* data, const long samples, const long max_lag, IOSTREAMD);

void write_result(double **array, long size, std::string outputfile_results_name, IOSTREAMD);
void write_multidim_result(double ***array, unsigned int dimens, long size, std::string outputfile_results_name, IOSTREAMD);
std::string bool2textMX(bool value);

void apply_baseline_correction(double* data, long samples);
