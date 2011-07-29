#include <cstring>
#include <cstdlib>
#include <iostream>
#include "te-datainit.h"

using namespace std;

// #define GSL_RANDOM_NUMBER_GENERATOR gsl_rng_default
#define GSL_RANDOM_NUMBER_GENERATOR gsl_rng_ranlxs2


int main()
{
  cout <<" ------ init ------ "<<endl;
  // initialize random number generator
  gsl_rng_env_setup();
  gsl_rng* GSLshuffler = gsl_rng_alloc(GSL_RANDOM_NUMBER_GENERATOR);
  gsl_rng_set(GSLshuffler, rand());
  
  
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
  long samples = 20*60*1000/tauImg;
  string model = "Leogang"; // "HowManyAreActive"; //"SpikeCount";
  bool OverrideRescalingQ = false;
  double std_noise = 0.03; //-1.;
  double fluorescence_saturation = 300; //-1.;
  double sigma_scatter = 0.15;
  double amplitude_scatter = 0.15;
  
  cout <<endl<<" ------ generating time series from spike data ------ "<<endl;
  double** xdata = generate_time_series_from_spike_data(inputfile_spiketimes,inputfile_spikeindices,size,\
    tauImg,samples,model,std_noise,fluorescence_saturation,-1.,50.,1000.,GSLshuffler);
  cout <<"-> "<<samples<<" samples generated."<<endl;
  // display_subset(xdata[0]);
  
  // cout <<endl<<" ------ guessing conditioning level ------ "<<endl;
  // double GlobalConditioningLevel = Magic_GuessConditioningLevel(xdata,size,samples);
  // cout <<"-> conditioning level is: "<<GlobalConditioningLevel<<endl;

  // cout <<endl<<" ------ simulating light scattering ------ "<<endl;
  // apply_light_scattering_to_time_series(xdata,size,samples,YAMLfilename,sigma_scatter,amplitude_scatter);
  // cout <<"-> done."<<endl;

  cout <<endl<<" ------ generating global signal ------ "<<endl;
  rawdata globalbins = 2;
  double GlobalConditioningLevel = 0.18;
  unsigned long* AvailableSamples = new unsigned long[globalbins];
  long StartSampleIndex = 0;
  long EndSampleIndex = samples-1;
  rawdata* xglobal = generate_discretized_global_time_series(xdata,size,samples,globalbins, \
    GlobalConditioningLevel,AvailableSamples,StartSampleIndex,EndSampleIndex);
  
  cout <<endl<<" ------ applying high-pass filter ------ "<<endl;
  apply_high_pass_filter_to_time_series(xdata,size,samples);
  // display_subset(xdata[0]);
  cout <<"-> done."<<endl;
  
  // cout <<endl<<" ------ guessing number of bins ------ "<<endl;
  // int nr_bins = Magic_GuessBinNumber(xdata,size,samples);
  // cout <<"-> number of bins is: "<<nr_bins<<endl;
  // 
  // cout <<endl<<" ------ discretizing time series ------ "<<endl;
  // rawdata** xdataHP = generate_discretized_version_of_time_series(xdata,size,samples,nr_bins);
  // display_subset(xdataHP[0]);

  /* cout <<endl<<" ------ test differential entropy (utility functions) ------ "<<endl;
  for(int i=1; i<=3; i++)
    cout <<"-> S_"<<i<<" = "<<SphericalUnitSurface(i)<<endl;
  gsl_vector* vec1 = gsl_vector_alloc(3);
  gsl_vector_set(vec1,0,0.5);
  gsl_vector_set(vec1,1,0.5);
  gsl_vector_set(vec1,2,0.5);
  gsl_vector* vec2 = gsl_vector_alloc(3);
  gsl_vector_set(vec2,0,0.5);
  gsl_vector_set(vec2,1,9.5);
  gsl_vector_set(vec2,2,0.5);
  cout <<"-> gsl_quicknorm(vec1,vec2,dim) = "<<gsl_quicknorm(vec1,vec2,3)<<endl;
  
  cout <<endl<<" ------ test differential entropy (0_Now and 1_Now) ------ "<<endl;
  gsl_vector** xtest = new gsl_vector*[samples];
  for (long s=0; s<samples; s++) {
    xtest[s] = gsl_vector_alloc(2);
    // fill vectors with 0_Now and 1_Now
    gsl_vector_set(xtest[s],0,xdata[0][s]);
    gsl_vector_set(xtest[s],1,xdata[1][s]);
  }
  double Hdiff = DifferentialEntropy(xtest,2,samples);
  cout <<"-> differential entropy is: "<<Hdiff<<endl;
  
  cout <<endl<<" ------ test differential entropy (0_Now and 0_Past) ------ "<<endl;
  gsl_vector_set(xtest[0],0,0.);
  gsl_vector_set(xtest[0],1,0.);
  for (long s=1; s<samples; s++) {
    // fill vectors with 0_Now and 0_Past
    gsl_vector_set(xtest[s],0,xdata[0][s]);
    gsl_vector_set(xtest[s],1,xdata[0][s-1]);
  }
  Hdiff = DifferentialEntropy(xtest,2,samples);
  cout <<"-> differential entropy is: "<<Hdiff<<endl; */
  
  // initialize random number generator
  gsl_rng_env_setup();
  gsl_rng* GSLrandom = gsl_rng_alloc(GSL_RANDOM_NUMBER_GENERATOR);
  gsl_rng_set(GSLrandom, rand());
  
  cout <<endl<<" ------ test geometric shuffling (part Ia) ------ "<<endl;
  long test_samples = 20;
  long* test_series = new long[test_samples];
  for(long t=0; t<test_samples; t++) test_series[t] = t;
  random_permutation(&test_series,test_samples,GSLrandom);
  for(long t=0; t<test_samples; t++) cout <<"test_series: (t,x) = ("<<t<<","<<test_series[t]<<")"<<endl;
  
  // cout <<endl<<" ------ test geometric shuffling (part Ib) ------ "<<endl;
  // // long test_samples = 20;
  // // long* test_series = new long[test_samples];
  // for(long t=0; t<test_samples; t++) test_series[t] = t;
  // geometric_permutation(test_series,test_samples,3);
  // for(long t=0; t<test_samples; t++) cout <<"test_series: (t,x) = ("<<t<<","<<test_series[t]<<")"<<endl;
  
  // cout <<endl<<" ------ test geometric shuffling (part II) ------ "<<endl;
  // int AutoCorrLength = 5;
  // long* shuffle = generate_random_geometric_permutation(samples,globalbins,xglobal,AutoCorrLength);
  // cout <<"s\tg\tsh"<<endl<<"----\t----\t----"<<endl;
  // for(long s=0; s<samples; s++) {
  //   cout <<s<<"\t"<<int(xglobal[s])<<"\t"<<shuffle[s]<<"\t"<<endl;
  // }

  cout <<endl<<" ------ test auto-correlation time scale ------ "<<endl;
  double tau = AutoCorrelationTimeScale(xdata[0],samples,100);
  cout <<"-> tau = "<<tau<<std::endl;
  
  cout <<endl<<" ------ deallocating memory ------ "<<endl;
  free_time_series_memory(xdata, size);
  gsl_rng_free(GSLrandom);  
  // free_time_series_memory(xdataHP, size);
  // for (long s=0; s<samples; s++) gsl_vector_free(xtest[s]);
  // delete[] xtest;  
  
  // cout <<endl<<" ------ testing YAML import ------ "<<endl;
  // read_positions_from_YAML("multi-topologies/Leogang/adjA_Leogang1.yaml", 100);
  
  cout <<"done."<<endl;
  return 0;
};
