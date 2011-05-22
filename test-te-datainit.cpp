#include <cstring>
#include <cstdlib>
#include <iostream>
#include "te-datainit.h"

using namespace std;


int main()
{
  cout <<" ------ init ------ "<<endl;
  
  // cout <<"testing minimum:"<<endl;
  // double test1[] = {0.2, 3.5, 4.4};
  // double test2[] = {1.2, 13.5, 14.4,-2.1};
  // double** test2d = new double*[2];
  // test2d[0] = &test1[0];
  // test2d[1] = &test2[0];
  // cout <<smallest(test2d, (unsigned int)2, (long)3)<<endl;
  
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

  // string inputfile_spikeindices = "multi-topologies/Leogang/s_index_10.dat";
  // string inputfile_spiketimes = "multi-topologies/Leogang/s_times_10.dat";
  string inputfile_spikeindices = "../Simulationen/NEST/calciumbursts2/Paris/LeogangTopology/s_index.dat";
  string inputfile_spiketimes = "../Simulationen/NEST/calciumbursts2/Paris/LeogangTopology/s_times.dat";
  string YAMLfilename = "../Simulationen/NEST/calciumbursts2/topology-Leogang.yml";
  cout <<endl<<" ------ loading files ------ "<<endl;
  cout <<"- \""<<inputfile_spikeindices<<"\""<<endl;
  cout <<"- \""<<inputfile_spiketimes<<"\""<<endl;
  unsigned int size = 100;
  unsigned int tauImg = 20;
  // long samples = 2*60*60*1000/tauImg;
  long samples = 10*1000/tauImg;
  string model = "Leogang"; // "HowManyAreActive"; //"SpikeCount";
  bool OverrideRescalingQ = false;
  double std_noise = 0.03; //-1.;
  double fluorescence_saturation = 300; //-1.;
  double sigma_scatter = 0.15;
  double amplitude_scatter = 0.15;
  
  cout <<endl<<" ------ generating time series from spike data ------ "<<endl;
  double** xdata = generate_time_series_from_spike_data(inputfile_spiketimes,inputfile_spikeindices,size,\
    tauImg,samples,model,std_noise,fluorescence_saturation,-1.,50.,1000.);
  cout <<"-> "<<samples<<" samples generated."<<endl;
  // display_subset(xdata[0]);
  
  cout <<endl<<" ------ guessing conditioning level ------ "<<endl;
  double GlobalConditioningLevel = Magic_GuessConditioningLevel(xdata,size,samples);
  cout <<"-> conditioning level is: "<<GlobalConditioningLevel<<endl;

  cout <<endl<<" ------ simulating light scattering ------ "<<endl;
  apply_light_scattering_to_time_series(xdata,size,samples,YAMLfilename,sigma_scatter,amplitude_scatter);
  cout <<"-> done."<<endl;

  cout <<endl<<" ------ guessing conditioning level (again) ------ "<<endl;
  GlobalConditioningLevel = Magic_GuessConditioningLevel(xdata,size,samples);
  cout <<"-> conditioning level is: "<<GlobalConditioningLevel<<endl;
  
  // cout <<endl<<" ------ guessing conditioning level (minima at zero test) ------ "<<endl;
  // Test_SetMinimaToZero(xdata,size,samples);
  // GlobalConditioningLevel = Magic_GuessConditioningLevel(xdata,size,samples);
  // cout <<"-> conditioning level is: "<<GlobalConditioningLevel<<endl;
  
  // cout <<endl<<" ------ applying high-pass filter ------ "<<endl;
  // apply_high_pass_filter_to_time_series(xdata,size,samples);
  // display_subset(xdata[0]);
  
  // cout <<endl<<" ------ guessing number of bins ------ "<<endl;
  // int nr_bins = Magic_GuessBinNumber(xdata,size,samples);
  // cout <<"-> number of bins is: "<<nr_bins<<endl;
  // 
  // cout <<endl<<" ------ discretizing time series ------ "<<endl;
  // rawdata** xdataHP = generate_discretized_version_of_time_series(xdata,size,samples,nr_bins);
  // display_subset(xdataHP[0]);
  
  cout <<endl<<" ------ deallocating memory ------ "<<endl;
  free_time_series_memory(xdata, size);
  // free_time_series_memory(xdataHP, size);
  
  // cout <<endl<<" ------ testing YAML import ------ "<<endl;
  // read_positions_from_YAML("multi-topologies/Leogang/adjA_Leogang1.yaml", 100);
  
  cout <<"done."<<endl;
  return 0;
};
