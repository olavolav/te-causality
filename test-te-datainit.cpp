#include <cstring>
#include <cstdlib>
#include <iostream>
#include "te-datainit.h"

#define SUBSET_LENGTH 100

using namespace std;

void display_subset(double* data);
void display_subset(rawdata* data);

int main()
{
  cout <<"start..."<<endl;
  
  // string inputfile_name = "../Simulationen/NEST/calciumbursts2/Paris/LeogangTopology/noisefree_rawCalcium_uchar_20ms.dat";
  // cout <<"loading file: '"<<inputfile_name<<"' ..."<<endl;
  // unsigned int size = 100;
  // long samples = 359951;
  // double input_scaling = 255./500.;
  // bool OverrideRescalingQ = false;
  // double std_noise = 10.; //-1.;
  // double fluorescence_saturation = -1.;
  // double cutoff = -1.;
  // double** xdata = load_time_series_from_binary_file(inputfile_name, size, samples, input_scaling, OverrideRescalingQ, std_noise, fluorescence_saturation, cutoff);
  // display_subset(xdata[0]);
  // apply_high_pass_filter_to_time_series(xdata, size, samples);
  // display_subset(xdata[0]);
  // unsigned int nr_bins = 10;
  // rawdata** xdataHP = generate_discretized_version_of_time_series(xdata, size, samples, nr_bins);
  // display_subset(xdataHP[0]);
  // 
  // cout <<"deallocating memory..."<<endl;
  // free_time_series_memory(xdata, size);
  // free_time_series_memory(xdataHP, size);

  string inputfile_spikeindices = "multi-topologies/Leogang/s_index_10.dat";
  string inputfile_spiketimes = "multi-topologies/Leogang/s_times_10.dat";
  cout <<"loading files:"<<endl;
  cout <<"- \""<<inputfile_spikeindices<<"\""<<endl;
  cout <<"- \""<<inputfile_spiketimes<<"\""<<endl;
  unsigned int size = 100;
  unsigned int tauImg = 20;
  long samples = SUBSET_LENGTH;
  string model = "Leogang"; // "HowManyAreActive"; //"SpikeCount";
  bool OverrideRescalingQ = false;
  double std_noise = 10.; //-1.;
  double fluorescence_saturation = -1.;
  double cutoff = -1.;
  double** xdata = generate_time_series_from_spike_data(inputfile_spiketimes, inputfile_spikeindices, size, tauImg, samples, model, std_noise, fluorescence_saturation, cutoff);
  cout <<"-> "<<samples<<" samples generated."<<endl;
  display_subset(xdata[0]);
  // apply_high_pass_filter_to_time_series(xdata, size, samples);
  // display_subset(xdata[0]);
  // unsigned int nr_bins = 10;
  // rawdata** xdataHP = generate_discretized_version_of_time_series(xdata, size, samples, nr_bins);
  // display_subset(xdataHP[0]);
  
  cout <<"deallocating memory..."<<endl;
  free_time_series_memory(xdata, size);
  // free_time_series_memory(xdataHP, size);
  
  cout <<"done."<<endl;
  return 0;
};

void display_subset(double* data)
{
  cout <<"displaying some subset of data points:"<<endl<<"{";
  for (long t=0; t<SUBSET_LENGTH; t++)
  {
    if (t>0) cout <<",";
    cout <<data[t];
  }
  cout <<"} (range "<<smallest(data,SUBSET_LENGTH)<<" – "<<largest(data,SUBSET_LENGTH)<<")"<<endl;
};
void display_subset(rawdata* data)
{
  cout <<"displaying some subset of data points:"<<endl<<"{";
  for (long t=0; t<SUBSET_LENGTH; t++)
  {
    if (t>0) cout <<",";
    cout <<int(data[t]);
  }
  cout <<"} (range "<<smallest(data,SUBSET_LENGTH)<<" – "<<largest(data,SUBSET_LENGTH)<<")"<<endl;
};
