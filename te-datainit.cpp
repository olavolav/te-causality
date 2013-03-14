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

#include "te-datainit.h"

#define FMODEL_SPIKECOUNT 1
#define FMODEL_HOWMANYAREACTIVE 2
#define FMODEL_LEOGANG 3
#define FMODEL_ERROR -1

#define HEIGHT_OF_ASCII_PLOTS 12

#define OUTPUTNUMBER_PRECISION 15

#undef SPIKE_INPUT_DATA_IS_BINARY
#undef TIME_SERIES_INPUT_DATA_IS_BINARY
#define NUMBER_OF_ROWS_TO_SKIP_IN_TIME_SERIES_INPUT_FILE 0 // 1 for Jordi's files
#define BASELINE_CORRECTION_BANDWIDTH 2000
#define SPEEDUP_BASELINE_CORRECTION

#define MAX_NUMBER_OF_BAD_DATA_LINES 20

// set output stream depending on wether SimKernel's sim.h is included
// (see also te-datainit.h)
#undef IOSTREAMH
#undef IOSTREAMC
#undef IOSTREAMV

#ifdef SIM_IO_H
  // SimKernel found.
  #define IOSTREAMH Sim& output
  #define IOSTREAMC output.io
  #define IOSTREAMV output
  #define IOSTREAMENDL Endl
#else
  // SimKernel not found, using standard output.
  #define IOSTREAMH std::ostream* output
  #define IOSTREAMC *output
  #define IOSTREAMV output
  #define IOSTREAMENDL std::endl
#endif

using namespace std;

double** load_time_series_from_file(std::string inputfile_name, unsigned int size, long samples, double input_scaling, bool OverrideRescalingQ, double std_noise, double fluorescence_saturation, double cutoff, gsl_rng* GSLrandom, IOSTREAMH)
{  
  // reserve and clear memory for result ("try&catch" is still missing!)
  double **xresult = NULL;
  char* in_from_file_array = NULL;
  double* tempdoublearray = NULL;
  try {
    xresult = new double*[size];
    for(unsigned int i=0; i<size; i++)
    {
      xresult[i] = NULL;
      xresult[i] = new double[samples];
      memset(xresult[i], 0, samples*sizeof(double));
    }
    in_from_file_array = new char[samples];
    tempdoublearray = new double[samples];
  }
  catch(...) {
    IOSTREAMC <<IOSTREAMENDL<<"error in load_time_series_from_file: cannot allocate enough memory!"<<IOSTREAMENDL;
    exit(1);
  }
  // open input file
  char* name = new char[inputfile_name.length()+1];
  strcpy(name,inputfile_name.c_str());
#ifdef TIME_SERIES_INPUT_DATA_IS_BINARY
  IOSTREAMC <<"-> setting up binary input ..."<<IOSTREAMENDL;
  std::ifstream inputfile(name, std::ios::binary);
#else
  IOSTREAMC <<"-> setting up plain text input ..."<<IOSTREAMENDL;
  std::ifstream inputfile(name);
#endif
  delete[] name;
  if (inputfile == NULL) {
    IOSTREAMC <<IOSTREAMENDL<<"error in load_time_series_from_file: cannot find input file!"<<IOSTREAMENDL;
    exit(1);
  }

  // test file length
#ifndef TIME_SERIES_INPUT_DATA_IS_BINARY
  long apparent_size = 1;
  long apparent_samples = 0;

  string line;
  int temp_pos;
  bool first_line = true;
  while (getline(inputfile, line)) {
    apparent_samples++;
    if(first_line) {
      while((temp_pos = line.find(",")) != std::string::npos) {
        apparent_size++;
        line.replace(temp_pos,1," ");
      }
    }
    first_line = false;
  }
  inputfile.clear();
  inputfile.seekg(0);
#if NUMBER_OF_ROWS_TO_SKIP_IN_TIME_SERIES_INPUT_FILE>0
  IOSTREAMC <<"Warning in load_time_series_from_file: Set to skip first ";
  IOSTREAMC <<NUMBER_OF_ROWS_TO_SKIP_IN_TIME_SERIES_INPUT_FILE<<" rows of each column."<<IOSTREAMENDL;
#endif
  apparent_size -= NUMBER_OF_ROWS_TO_SKIP_IN_TIME_SERIES_INPUT_FILE;    // because 1st line is sample number
  
  IOSTREAMC <<"-> it appears that the file contains "<<apparent_size<<" nodes and "<<apparent_samples;
  IOSTREAMC <<" samples each."<<IOSTREAMENDL;
  
  if(apparent_size < 1) {
    IOSTREAMC <<IOSTREAMENDL<<"error in load_time_series_from_file: could not detect number of nodes in file!"<<IOSTREAMENDL;
    inputfile.close();
    exit(1);
  }
  if(apparent_samples < 2) {
    IOSTREAMC <<IOSTREAMENDL<<"error in load_time_series_from_file: could not detect number of samples in file!"<<IOSTREAMENDL;
    inputfile.close();
    exit(1);
  }
  if(apparent_size != size) {
    IOSTREAMC <<IOSTREAMENDL<<"error in load_time_series_from_file: number of nodes in file does not match given size!"<<IOSTREAMENDL;
    inputfile.close();
    exit(1);
  }
  if(apparent_samples < samples) {
    IOSTREAMC <<IOSTREAMENDL<<"error in load_time_series_from_file: number of lines in file is lower than given sample number!"<<IOSTREAMENDL;
    inputfile.close();
    exit(1);
  }
#else
  inputfile.seekg(0,std::ios::end);
  if(long(inputfile.tellg()) != size*samples)
  {
    IOSTREAMC <<IOSTREAMENDL<<"error in load_time_series_from_file: file length of input does not match given parameters!"<<IOSTREAMENDL;
    exit(1);
  }
  inputfile.seekg(0,std::ios::beg);
#endif

  // import and rescale data  
  IOSTREAMC <<"reading data..."<<IOSTREAMENDL;
#ifdef TIME_SERIES_INPUT_DATA_IS_BINARY
  for(int j=0; j<size; j++)
  {
    inputfile.read(in_from_file_array, samples);

    // OverrideRescalingQ = true
    // Dies ignoriert also "appliedscaling", "noise", "HighPassFilterQ" und "cutoff"
    // Therefore, "bins" takes the role of an upper cutoff
    if (OverrideRescalingQ)
      for(long k=0; k<samples; k++)
        // xdata[j][k] = in_from_file_array[k];
        xresult[j][k] = double(in_from_file_array[k]);
    else     // OverrideRescalingQ = false
    {
      for (long k=0; k<samples; k++)
      {
        // transform to unsigned notation
        tempdoublearray[k] = double(in_from_file_array[k]);
        if (in_from_file_array[k]<0) tempdoublearray[k] += 256.;
        
        // transform back to original signal (same as in Granger case)
        tempdoublearray[k] /= input_scaling;
        // assuming a saturation with hill function of order 1
        if (fluorescence_saturation > 0.)
          tempdoublearray[k] = tempdoublearray[k]/(tempdoublearray[k] + fluorescence_saturation);
        // adding noise
        if (std_noise > 0.)
          tempdoublearray[k] += gsl_ran_gaussian(GSLrandom,std_noise);          
        // apply cutoff
        if ((cutoff>0)&&(tempdoublearray[k]>cutoff)) tempdoublearray[k] = cutoff;
      }
      
    }
    memcpy(xresult[j],tempdoublearray,samples*sizeof(double));
  }
#else
  int next_pos;
  int length;

  for (long tt=0; tt<samples; tt++) {
    // if((tt%5000)==0) {
    //   cout <<"debug: reading sample #"<<tt<<" ..."<<endl;
    // }
    getline(inputfile, line);
    length = line.length();
    // cout <<"debug: read = "<<line<<endl;
    
    // before accepting the values, see if the number of commas matches
    temp_pos = -1;
    int comma_count = 0;
    while((temp_pos = line.find(",",temp_pos+1)) != std::string::npos) {
      comma_count++;
    }
    // for (int i=0; i<size+NUMBER_OF_ROWS_TO_SKIP_IN_TIME_SERIES_INPUT_FILE+1; i++) {
    //   temp_pos = line.find(",",temp_pos+1);
    //   if(temp_pos!=std::string::npos) comma_count++;
    // }
    long skipped_rows_count = 0;
    if(comma_count+1 != size+NUMBER_OF_ROWS_TO_SKIP_IN_TIME_SERIES_INPUT_FILE) {
      IOSTREAMC <<"warning in load_time_series_from_file: skipping line #"<<tt;
      IOSTREAMC <<", because the number of entries is "<<comma_count<<" instead of ";
      IOSTREAMC <<size+NUMBER_OF_ROWS_TO_SKIP_IN_TIME_SERIES_INPUT_FILE<<"!"<<IOSTREAMENDL;
      
      // copying old values
      assert(tt>0);
      for (int i=0; i<size; i++) {
        xresult[i][tt] = xresult[i][tt-1];
      }
      skipped_rows_count++;
      assert(skipped_rows_count<MAX_NUMBER_OF_BAD_DATA_LINES);
    }
    else {
      temp_pos = -1;
      for (int i=0; i<size+NUMBER_OF_ROWS_TO_SKIP_IN_TIME_SERIES_INPUT_FILE; i++) {
        next_pos = line.find(",",temp_pos+1);
        if(next_pos==std::string::npos) {
          next_pos = length - 1;
        }
        // line.replace(next_pos,1," ");
        if (i >= NUMBER_OF_ROWS_TO_SKIP_IN_TIME_SERIES_INPUT_FILE) {
          xresult[i-NUMBER_OF_ROWS_TO_SKIP_IN_TIME_SERIES_INPUT_FILE][tt] = \
            atof(line.substr(temp_pos+1,next_pos-temp_pos+1).c_str());
          // cout <<"debug: xresult[i][tt] = "<<xresult[i-NUMBER_OF_ROWS_TO_SKIP_IN_TIME_SERIES_INPUT_FILE][tt]<<endl;
        }
        temp_pos = next_pos;
      }
      // exit(1);
      // if((tt%50000)==0) {
      //   cout <<"debug: sample #"<<tt<<": data of first three is: "<<xresult[0][tt]<<", "<<xresult[1][tt]<<", "<<xresult[2][tt]<<endl;
      // }
    }
  }
#endif
  inputfile.close();
  
  // apply a moning window correction (calculation see 24.10.11)
  // const long mwa_sigma = 2;
  // const long athird = samples/mwa_sigma; // implicit floor
  // IOSTREAMC <<"HACK WARNING: Moving window averaging activated, with a width of sigma = "<<mwa_sigma<<IOSTREAMENDL;
  // long shift, first, last, s2;
  // double* copy_array = NULL;
  // copy_array = new double[samples];
  // for (int i=0; i<size; i++) {
  //   memcpy(copy_array,xresult[i],samples*sizeof(double));
  //   for (long s=0; s<samples; s++) {
  //     shift = s/athird; // implicit floor
  //     if (shift < mwa_sigma) {
  //       s2 = s - shift*athird;
  //       first = mwa_sigma*s2 + shift;
  //       last = min(samples-1, mwa_sigma*(s2+1) + shift -1);
  //       // if (i==0) cout <<"debug: s="<<s<<": first="<<first<<", last="<<last<<endl;
  //       xresult[i][s] = mean(copy_array, first, last);
  //     }
  //   }
  // }
  // delete[] copy_array;

  // baseline correction
  IOSTREAMC <<"applying baseline correction..."<<IOSTREAMENDL;
  for (int i=0; i<size; i++) {
    apply_baseline_correction(xresult[i],samples);
  }
  // std::cout <<"debug: result for first node: ";
  // for (long tt=0; tt<samples; tt+=(long)(samples/13.)) {
  //   std::cout <<xresult[0][tt]<<std::endl;
  // }
  // exit(0);
  
  // applying offset such that averaged signal is never negative
  IOSTREAMC <<"applying offset such that averaged signal is never negative..."<<IOSTREAMENDL;
  double min_xmean = std::numeric_limits<double>::max();
  double temp_xmean;
  // first, find lowest xmean value...
  for(long tt=0; tt<samples; tt++) {
    temp_xmean = 0.;
    for (int i=0; i<size; i++) {
      temp_xmean += xresult[i][tt];
    }
    temp_xmean /= size;
    if(temp_xmean < min_xmean) min_xmean = temp_xmean;
  }
  // then, apply resulting offset
  for(long tt=0; tt<samples; tt++) {
    for (int i=0; i<size; i++) {
      xresult[i][tt] -= min_xmean;
    }
  }

  // // determine available samples per globalbin for TE normalization later
  // memset(AvailableSamples, 0, globalbins*sizeof(unsigned long));
  // for (unsigned long t=StartSampleIndex; t<EndSampleIndex; t++)
  //  AvailableSamples[xglobal[t]]++;
  // 
  // if (EqualSampleNumberQ || (MaxSampleNumberPerBin>0))
  // {
  //  unsigned long maxsamples = ULONG_MAX;
  //  for (rawdata g=0; g<globalbins; g++)
  //    if (AvailableSamples[g]<maxsamples) maxsamples = AvailableSamples[g];
  //  IOSTREAMC <<"DEBUG: maxsamples = "<<maxsamples<<IOSTREAMENDL;
  //  
  //  if ((MaxSampleNumberPerBin>maxsamples) && !EqualSampleNumberQ)
  //    maxsamples = MaxSampleNumberPerBin;
  //  if ((MaxSampleNumberPerBin<maxsamples)&&(MaxSampleNumberPerBin>0))
  //    maxsamples = MaxSampleNumberPerBin;
  //  IOSTREAMC <<"DEBUG: cut to maxsamples = "<<maxsamples<<IOSTREAMENDL;
  //  
  //  unsigned long* AlreadySelectedSamples = new unsigned long[globalbins];
  //  memset(AlreadySelectedSamples, 0, globalbins*sizeof(unsigned long));
  //  for (unsigned long t=StartSampleIndex; t<EndSampleIndex; t++)
  //    if ((++AlreadySelectedSamples[xglobal[t]])>maxsamples)
  //      xglobal[t] = globalbins; // ..and therefore exclude from calculation
  // 
  //  // re-determine available samples per globalbin (inefficient)
  //  memset(AvailableSamples, 0, globalbins*sizeof(unsigned long));
  //  for (unsigned long t=StartSampleIndex; t<EndSampleIndex; t++)
  //    if (xglobal[t]<globalbins) AvailableSamples[xglobal[t]]++;
  //    
  //  delete[] AlreadySelectedSamples;
  // }
  
  // free allocated memory
  // gsl_rng_free(GSLrandom);    
  delete[] in_from_file_array;
  // delete[] xglobaltemp;
  delete[] tempdoublearray;
  // delete[] tempdoublearraycopy;
  
  return xresult;
};



rawdata* generate_discretized_global_time_series(double** time_series, unsigned int size, long samples, unsigned int globalbins, double GlobalConditioningLevel, unsigned long* AvailableSamples, long StartSampleIndex, long EndSampleIndex, bool EqualSampleNumberQ, long MaxSampleNumberPerBin, IOSTREAMH)
{
  rawdata* xglobal = new rawdata[samples];
  memset(xglobal, 0, samples*sizeof(rawdata));
  double* xglobaltemp = generate_mean_time_series(time_series, size, samples);
    
  // EVIL SAALBACH HACK FOR TIME CODE GLOBAL SIGNAL: -------------------------------------------- !!!!!!!!!!
  // xglobaltemp[0] = 0.;
  // for (unsigned long t=0; t<samples; t++)
  //  xglobaltemp[t] = double(int(t)%int(60*24/tauF));
  //  // xglobaltemp[t] = double(mod(t,60*24/tauF));
  // IOSTREAMC <<"DEBUG OF EVIL TIME CODE HACK: last globaltemp value = "<<xglobaltemp[samples-1]<<IOSTREAMENDL;
  
  if (GlobalConditioningLevel > 0.)
  {
    unsigned long below = 0;
    for (long t=0; t<samples; t++)
    {
      if (xglobaltemp[t] > GlobalConditioningLevel) xglobal[t] = 1;
      else
      {
        xglobal[t] = 0;
        below++;
      }
    }
    IOSTREAMC <<" -> global conditioning level "<<GlobalConditioningLevel<<": "<<(100.*below)/samples;
    IOSTREAMC <<"% are below threshold. "<<IOSTREAMENDL;
  }
  else discretize(xglobaltemp,xglobal,samples,globalbins);
  
  // // EVIL HACK FOR SHIFTED GLOBAL BINS: -------------------------------------------- !!!!!!!!!!
  // IOSTREAMC <<"Warning: Evil global bin shifting hack enabled!!"<<IOSTREAMENDL;
  // const double gminlevel = -0.25;
  // const double gmaxlevel = 0.85;
  // discretize(xglobaltemp,xglobal,gminlevel,gmaxlevel,samples,globalbins);

  // EVIL HACK FOR SHIFTED GLOBAL BINS v2: -------------------------------------------- !!!!!!!!!!
  // (calculation see 03.10.12)
  // IOSTREAMC <<"Warning: Evil global bin shifting hack v2 enabled!!"<<IOSTREAMENDL;
  // const double gmaxlevel = largest(xglobaltemp, samples);
  // IOSTREAMC <<" -> Found maximum of xmean of x = "<<gmaxlevel<<IOSTREAMENDL;
  // const double ghistopeak = Util_FindPeakInHistogram(xglobaltemp, samples, smallest(xglobaltemp, samples), gmaxlevel, 60);
  // IOSTREAMC <<" -> Found peak in histogram at x = "<<ghistopeak<<IOSTREAMENDL;
  // const double gminlevel = (double(globalbins)*ghistopeak - gmaxlevel)/(double(globalbins)-1.0);
  // IOSTREAMC <<" -> Set new lower bound of histogram to x = "<<gminlevel<<IOSTREAMENDL;
  // discretize(xglobaltemp,xglobal,gminlevel,gmaxlevel,samples,globalbins);
  
  // determine available samples per globalbin for TE normalization later
  memset(AvailableSamples, 0, globalbins*sizeof(unsigned long));
  for (unsigned long t=StartSampleIndex; t<=EndSampleIndex; t++)
    AvailableSamples[xglobal[t]]++;

  if (EqualSampleNumberQ || (MaxSampleNumberPerBin>0)) {
    IOSTREAMC <<"Warning: Sample number overrides enabled!"<<IOSTREAMENDL;
    long maxsamples = LONG_MAX;
    for (rawdata g=0; g<globalbins; g++)
      if (AvailableSamples[g]<maxsamples) maxsamples = AvailableSamples[g];
    IOSTREAMC <<"DEBUG: maxsamples = "<<maxsamples<<IOSTREAMENDL;

    if ((MaxSampleNumberPerBin>maxsamples) && !EqualSampleNumberQ)
      maxsamples = MaxSampleNumberPerBin;
    if ((MaxSampleNumberPerBin<maxsamples)&&(MaxSampleNumberPerBin>0))
      maxsamples = MaxSampleNumberPerBin;
    IOSTREAMC <<"DEBUG: Cut to maxsamples = "<<maxsamples<<IOSTREAMENDL;

    unsigned long* AlreadySelectedSamples = new unsigned long[globalbins];
    memset(AlreadySelectedSamples, 0, globalbins*sizeof(unsigned long));
    for (unsigned long t=StartSampleIndex; t<EndSampleIndex; t++)
      if ((++AlreadySelectedSamples[xglobal[t]])>maxsamples)
        xglobal[t] = globalbins; // ..and therefore exclude from calculation

    // re-determine available samples per globalbin (inefficient)
    memset(AvailableSamples, 0, globalbins*sizeof(unsigned long));
    for (unsigned long t=StartSampleIndex; t<EndSampleIndex; t++)
      if (xglobal[t]<globalbins) AvailableSamples[xglobal[t]]++;
 
    delete[] AlreadySelectedSamples;
  }

  free_time_series_memory(xglobaltemp);
  return xglobal;
}

rawdata** generate_discretized_version_of_time_series(double** in, unsigned int size, long nr_samples, unsigned int nr_bins)
{
  rawdata** xout;
  xout = new rawdata*[size];
  for(unsigned int ii=0; ii<size; ii++)
  {
    xout[ii] = new rawdata[nr_samples];
    discretize(in[ii], xout[ii], nr_samples, nr_bins);
  }  
  return xout;
};

void discretize(double* in, rawdata* out, long nr_samples, unsigned int nr_bins)
{
  discretize(in,out,smallest(in,nr_samples),largest(in,nr_samples),nr_samples,nr_bins);
};
void discretize(double* in, rawdata* out, double min, double max, long nr_samples, unsigned int nr_bins)
{
  double xstepsize = (max-min)/nr_bins;

  int xint;
  for (unsigned long t=0; t<nr_samples; t++)
    out[t] = discretize(in[t],min,max,nr_bins);
};
rawdata discretize(double in, double min, double max, unsigned int nr_bins)
{
  assert(max>min);
  assert(nr_bins>0);
  // correct discretization according to 'te-test.nb'
  // incorporated later: double xstepsize = (max-min)/nr_bins;
  // IOSTREAMC <<"max = "<<max<<IOSTREAMENDL;
  // IOSTREAMC <<"min = "<<min<<IOSTREAMENDL;
  // IOSTREAMC <<"bins here = "<<nr_bins<<IOSTREAMENDL;
  // IOSTREAMC <<"stepsize = "<<xstepsize<<IOSTREAMENDL;

  int xint;

  // assert(in<=max); ...does not have to be true, and does not matter, data is included in highest bin then
  // assert(in>=min);
  if (in>=max) xint = nr_bins-1;
  else
  {
    if (in<=min) xint = 0;
    // with stepsize variable: else xint = (int)((in-min)/xstepsize);
    // without:
    else xint = (int)((in-min)*double(nr_bins)/(max-min));
  }
  if (xint >= nr_bins) xint = nr_bins-1; // need to have this for some silly numerical reason...

  assert((xint>=0)&&(rawdata(xint)<nr_bins)); // ...just to be sure...
  return rawdata(xint);
};

// Orlandi: Adding option for predefined binning limits
rawdata** generate_discretized_version_of_time_series(double** in, unsigned int size, long nr_samples, std::vector<double> binEdges)
{
  rawdata** xout;
  xout = new rawdata*[size];
  for(unsigned int ii=0; ii<size; ii++)
  {
    xout[ii] = new rawdata[nr_samples];
    discretize(in[ii], xout[ii], nr_samples, binEdges);
  }  
  return xout;
};

void discretize(double* in, rawdata* out, long nr_samples, std::vector<double> binEdges)
{
  int xint;
  for (unsigned long t=0; t<nr_samples; t++)
    out[t] = discretize(in[t], binEdges);
};

// For now the bottom and top edges are doing nothing, they act like (-inf and inf)
rawdata discretize(double in, std::vector<double> binEdges)
{
  // By default set it to the top bin
  int xint = binEdges.size()-2;

  // Correct to the right bin
  for(std::vector<double>::size_type i = 1; i != binEdges.size(); i++)
  {
    if(in < binEdges[i])
    {
      xint = i-1;
      break;
    }
  }
  assert((xint>=0)&&(rawdata(xint)<(binEdges.size()-1))); // ...just to be sure...*/
  return rawdata(xint);
};

void apply_high_pass_filter_to_time_series(double** time_series, unsigned int size, long nr_samples)
{
  for(unsigned int ii=0; ii<size; ii++)
    apply_high_pass_filter_to_time_series(time_series[ii], nr_samples);
};
void apply_high_pass_filter_to_time_series(double* time_series, long nr_samples)
{
  double* arraycopy = new double[nr_samples];

  // of course, this is just a difference signal, so not really filtered
  memcpy(arraycopy,time_series,nr_samples*sizeof(double));
  time_series[0] = 0.;
  for(long k=1; k<nr_samples; k++)
    time_series[k] = arraycopy[k] - arraycopy[k-1];

  delete[] arraycopy;
};

double** generate_time_series_from_spike_data(std::string inputfile_spiketimes, std::string inputfile_spikeindices, unsigned int size, unsigned int tauImg, long samples, std::string fluorescence_model, double std_noise, double fluorescence_saturation, double cutoff, double DeltaCalciumOnAP, double tauCa, gsl_rng* GSLrandom, IOSTREAMH)
{
  // reserve and clear memory for result ("try&catch" is still missing!)
  double **xresult = new double*[size];
  for(unsigned int i=0; i<size; i++)
  {
    xresult[i] = new double[samples];
    memset(xresult[i], 0, samples*sizeof(double));
  }
  
  // open files
  char* nameI = new char[inputfile_spikeindices.length()+1];
  strcpy(nameI,inputfile_spikeindices.c_str());
#ifdef SPIKE_INPUT_DATA_IS_BINARY
  std::ifstream inputfileI(nameI, std::ios::binary);
#else
  std::ifstream inputfileI(nameI);
#endif
  if (inputfileI == NULL) {
    IOSTREAMC <<IOSTREAMENDL<<"error: cannot find spike indices file!"<<IOSTREAMENDL;
    exit(1);
  }
  delete[] nameI;
  
  char* nameT = new char[inputfile_spiketimes.length()+1];
  strcpy(nameT,inputfile_spiketimes.c_str());
#ifdef SPIKE_INPUT_DATA_IS_BINARY
  IOSTREAMC <<"-> setting up binary input"<<IOSTREAMENDL;
  std::ifstream inputfileT(nameT, std::ios::binary);
#else
  IOSTREAMC <<"-> setting up plain text input"<<IOSTREAMENDL;
  std::ifstream inputfileT(nameT);
#endif
  if (inputfileT == NULL) {
    IOSTREAMC <<IOSTREAMENDL<<"error: cannot find spike times file!"<<IOSTREAMENDL;
    exit(1);
  }
  delete[] nameT;
  
  // determine file length, then allocate memory
  long nr_spikes = 0;
#ifdef SPIKE_INPUT_DATA_IS_BINARY
  inputfileI.seekg(0,std::ios::end);
  nr_spikes = inputfileI.tellg()/sizeof(int);
  inputfileI.seekg(0,std::ios::beg);
#else
  string line;
  long tempsize = 0;
  // while (!inputfileI.eof()) {
  while (getline(inputfileI, line)) {
    // cout <<"not the end of file!"<<flush;
    // getline(inputfileI, line);
    // cout <<"readin:"<<line<<endl;
    tempsize++;
  }
  nr_spikes = tempsize;
  inputfileI.clear();
  inputfileI.seekg(0);
#endif
  IOSTREAMC <<"-> number of spikes in index file: "<<nr_spikes<<IOSTREAMENDL;
  int* xindex = new int[nr_spikes];
  double* xtimes = new double[nr_spikes];
  
  // read spike data
#ifdef SPIKE_INPUT_DATA_IS_BINARY
  inputfileI.read((char*)xindex, nr_spikes*sizeof(int));
  inputfileT.read(reinterpret_cast<char*>(xtimes), nr_spikes*sizeof(double));
#else
  for (long tt=0; tt<nr_spikes; tt++) {
    getline(inputfileI, line);
    // cout <<"read="<<line<<endl;
    xindex[tt] = atoi(line.c_str());
    getline(inputfileT, line);
    // cout <<"read="<<line<<endl;
    xtimes[tt] = atof(line.c_str());
  }
#endif
  
  // close files
  inputfileI.close();
  inputfileT.close();

  // debug output
  // for(long t=0; t<min((long)20,nr_spikes); t++)
  //   IOSTREAMC <<"DEBUG: xindex = "<<xindex[t]<<", xtimes = "<<xtimes[t]<<IOSTREAMENDL;
  // exit(0);

  // test if read data appears valid
  for (long tt=0; tt<nr_spikes; tt++) {
    assert((xindex[tt]>=0)&&(xindex[tt]<size)); // indices are in allowed range
    if(tt>0) assert(xtimes[tt]>=xtimes[tt-1]); // spike times are an ordered sequence
    if(tt<nr_spikes-1) assert(xtimes[tt]<=xtimes[tt+1]);
  }
  
  // choose switch key for the fluorescence model
  int fluorescence_model_key = FMODEL_ERROR;
  if (fluorescence_model == "SpikeCount") fluorescence_model_key = FMODEL_SPIKECOUNT;
  if (fluorescence_model == "HowManyAreActive") fluorescence_model_key = FMODEL_HOWMANYAREACTIVE;
  if (fluorescence_model == "Leogang") fluorescence_model_key = FMODEL_LEOGANG;
  if(fluorescence_model_key == FMODEL_ERROR) {
    IOSTREAMC <<IOSTREAMENDL<<"error: unknown fluorescence model!"<<IOSTREAMENDL;
    exit(1);
  }
  
  // generate fluorescence data
  IOSTREAMC <<"-> generate fluorescence data using model: "<<fluorescence_model<<IOSTREAMENDL;
  // const int int_tauF = (int)round(tauF); // in ms
  unsigned long startindex = 1;
  unsigned long endindex = 0; // therefore, we miss the first spike!
  long dataindex = 0;
  unsigned long tinybit_spikenumber;
  // unsigned long ttExactMS = 0;
  for (unsigned long ttExactMS=0; ttExactMS<tauImg*samples; ttExactMS+=tauImg)
  {
    // determine starting and ending spike index of current frame
    while ((endindex+1<nr_spikes)&&(xtimes[endindex+1]<=ttExactMS+tauImg))
      endindex++;
    tinybit_spikenumber = std::max(endindex-startindex+1,(unsigned long)0);
    
    // IOSTREAMC <<"DEBUG: ttExactMS = "<<ttExactMS<<", startindex = "<<startindex<< \
    //   ", endindex = "<<endindex<<", tinybit_spikenumber = "<<tinybit_spikenumber<<IOSTREAMENDL;

    for (int ii=0; ii<size; ii++)
    {
      switch (fluorescence_model_key)
      {
        case FMODEL_SPIKECOUNT:
          if(tinybit_spikenumber>0)
            xresult[ii][dataindex] = double(count(xindex,startindex,endindex,ii));
            // test: xresult[ii][dataindex] = double(tinybit_spikenumber);
          break;
          
        case FMODEL_HOWMANYAREACTIVE:
          if(tinybit_spikenumber>0)
            xresult[ii][dataindex] = double(has_index(xindex,startindex,endindex,ii));
          break;
          
        case FMODEL_LEOGANG:
          xresult[ii][dataindex] = (1.-double(tauImg)/tauCa)*xresult[ii][std::max(dataindex-1,long(0))] + \
            DeltaCalciumOnAP*double(count(xindex,startindex,endindex,ii));
          break;
          
        // default:
        //   IOSTREAMC <<"error in generate_time_series_from_spike_data: invalid fluorescence model"<<IOSTREAMENDL;
        //   exit(1);
      }
    }
    if(startindex <= endindex)
      startindex = 1 + endindex;
    dataindex++;
  }

  // apply saturation (Hill function of order 1 as usual)
  if(fluorescence_saturation > 0.)
    for (unsigned int ii=0; ii<size; ii++)
      for (long tt=0; tt<samples; tt++)
        xresult[ii][tt] = xresult[ii][tt]/(xresult[ii][tt]+fluorescence_saturation);

  // apply additive noise term
  if(std_noise > 0.)
  {
    // initialize random number generator
    // gsl_rng* GSLrandom;
    // gsl_rng_env_setup();
    // GSLrandom = gsl_rng_alloc(GSL_RANDOM_NUMBER_GENERATOR);
    // gsl_rng_set(GSLrandom, 1234);
    
    for (unsigned int ii=0; ii<size; ii++)
      for (long tt=0; tt<samples; tt++)
        xresult[ii][tt] += gsl_ran_gaussian(GSLrandom,std_noise);
        
    // free allocated memory
    // gsl_rng_free(GSLrandom);    
  }
  
  delete[] xindex;
  delete[] xtimes;
  
  return xresult;
};

unsigned long count(int* array, unsigned long starti, unsigned long endi, int what)
{
  unsigned long occur = 0;
  for (unsigned long i=starti; i<=endi; i++)
    if (array[i] == what) occur++;
  return occur;
};
bool has_index(int* array, unsigned long starti, unsigned long endi, int what)
{
  for (unsigned long i=starti; i<=endi; i++)
    if (array[i] == what) return true;
  return false;
};

double smallest(double* array, const long length)
{
  double min = array[0];
  for (long i=1; i<length; i++)
    if(array[i]<min) min = array[i];

  return min;
};
double largest(double* array, const long length)
{
  double max = array[0];
  for (long i=1; i<length; i++)
    if(array[i]>max) max = array[i];

  return max;
};
rawdata smallest(rawdata* array, const long length)
{
  rawdata min = array[0];
  for (long i=1; i<length; i++)
    if(array[i]<min) min = array[i];

  return min;
};
rawdata largest(rawdata* array, const long length)
{
  rawdata max = array[0];
  for (long i=1; i<length; i++)
    if(array[i]>max) max = array[i];

  return max;
};

double smallest(double** array, const unsigned int size, const long length)
{
  double min = array[0][0];
  for (unsigned int i=1; i<size; i++)
    min = std::min(min,smallest(array[i],length));

  return min;
};
double largest(double** array, const unsigned int size, const long length)
{
  double max = array[0][0];
  for (unsigned int i=1; i<size; i++)
    max = std::max(max,largest(array[i],length));

  return max;
};

double total(double* array, const long length)
{
  return total(array,0,length-1);
};
double total(double* array, const long first, const long last)
{
  if (last < first) return 0.;
  if (first == last) return array[first];
  
  double sum = 0.;
  for (long i=first; i<=last; i++)
    sum += array[i];

  return sum;
};

double mean(double* array, const long first, const long last)
{
  if (last < first) return 0.;
  if (first == last) return array[first];
  
  return total(array,first,last)/double(last-first+1);
}
double mean(double* array, const long length) {
  return mean(array,0,length-1);
};

double variance(double* array, const long first, const long last) {
  if (last <= first) return 0.0;
  double mu = mean(array, first, last);
  double sum_of_squares = 0.0;
  
  for(long i=first; i<=last; i++) {
    sum_of_squares += (array[i] - mu) * (array[i] - mu);
  }
  
  return sum_of_squares/double(last-first+1);
}
double variance(double* array, const long length) {
  return variance(array,0,length-1);
}

double standard_deviation(double* array, const long first, const long last) {
  return sqrt( variance(array, first, last) );
}
double standard_deviation(double* array, const long length) {
  return standard_deviation(array,0,length-1);
}

double* generate_mean_time_series(double** data, unsigned int size, long samples)
{
  double* xglobaltemp = new double[samples];
  memset(xglobaltemp, 0, samples*sizeof(double));

  for (long t=0; t<samples; t++)
  {
    for (unsigned int ii=0; ii<size; ii++)
      xglobaltemp[t] += data[ii][t];
    xglobaltemp[t] /= size;
  }
  
  return xglobaltemp;
}

void free_time_series_memory(double** xresult, unsigned int size)
{
  for(unsigned int ii=0; ii<size; ii++)
    free_time_series_memory(xresult[ii]);
  delete[] xresult;
};
void free_time_series_memory(double* xresult)
{
  delete[] xresult;
};
void free_time_series_memory(rawdata** xresult, unsigned int size)
{
  for(unsigned int ii=0; ii<size; ii++)
    free_time_series_memory(xresult[ii]);
  delete[] xresult;
};
void free_time_series_memory(rawdata* xresult)
{
  delete[] xresult;
};

void display_subset(double* data, int length, IOSTREAMH)
{
  // IOSTREAMC <<"displaying some subset of data points:"<<IOSTREAMENDL;
  IOSTREAMC <<"{";
  for (long t=0; t<length; t++)
  {
    if (t>0) IOSTREAMC <<",";
    IOSTREAMC <<data[t];
  }
  IOSTREAMC <<"} (range "<<smallest(data,length)<<" – "<<largest(data,length)<<")"<<IOSTREAMENDL;
};
void display_subset(rawdata* data, int length, IOSTREAMH)
{
  // IOSTREAMC <<"displaying some subset of data points:"<<IOSTREAMENDL;
  IOSTREAMC <<"{";
  for (long t=0; t<length; t++)
  {
    if (t>0) IOSTREAMC <<",";
    IOSTREAMC <<int(data[t]);
  }
  IOSTREAMC <<"} (range "<<int(smallest(data,length))<<" – "<<int(largest(data,length))<<")"<<IOSTREAMENDL;
};

int Magic_GuessBinNumber(double** data, const unsigned int size, const long samples)
{
  double range, std;  
  double meanbins = 0.;
  for(unsigned int i=0; i<size; i++)
  {
    range = largest(data[i],samples)-smallest(data[i],samples);
    assert(range > 0.);
    std = sqrt(gsl_stats_variance(data[i],1,samples));
    meanbins += 1*range/std; // old code: 2*std/range
  }
  meanbins /= size;
  cout <<"debug: meanbins = "<<meanbins<<endl;
  
  return std::max(2,int(round(meanbins)));
  // return std::max(2,int(meanbins));
};

double Magic_GuessConditioningLevel(double** data, const unsigned int size, const long samples, IOSTREAMH)
{
  int histo_bins = int(std::max(4.,std::min(200.,sqrt(samples))));
  IOSTREAMC <<" -> number of bins for histogram: "<<histo_bins<<IOSTREAMENDL;
  double xmeanmin, xmeanmax;
  double xresultlevel = -1.;
  
  double* xmean = generate_mean_time_series(data,size,samples);
  xmeanmin = smallest(xmean,samples);
  xmeanmax = largest(xmean,samples);
  // IOSTREAMC <<"-> xmeanmin = "<<xmeanmin<<", xmeanmax = "<<xmeanmax<<IOSTREAMENDL;
  // IOSTREAMC <<" -> beginning of <f>: "<<IOSTREAMENDL;
  // display_subset(xmean,5,IOSTREAMV);
  
  // find maximum, which we assume comes from the noise peak
  double x_ymax = Util_FindPeakInHistogram(xmean,samples,xmeanmin,xmeanmax,histo_bins);
  IOSTREAMC <<" -> identified peak at <f> = "<<x_ymax<<IOSTREAMENDL;
  
  // for(int i=0; i<histo_bins; i++)
  //   IOSTREAMC <<"histo <f> = "<<xmeanmin+(double(i)+0.5)*(xmeanmax-xmeanmin)/double(histo_bins)<<": count = "<<histo[i]<<IOSTREAMENDL;
  
  // PlotHistogramInASCII(xmean,samples,xmeanmin,x_ymax+0.1,"<f>","#(<f>)",IOSTREAMV);
  IOSTREAMC <<IOSTREAMENDL<<" -> log histogram of complete range of <f>:"<<IOSTREAMENDL;
  PlotLogHistogramInASCII(xmean,samples,xmeanmin,xmeanmax,"<f>","log #(<f>)",IOSTREAMV);
  // IOSTREAMC <<IOSTREAMENDL<<" -> log histogram of the right tail of <f>:"<<IOSTREAMENDL;
  // PlotLogHistogramInASCII(xmean,samples,x_ymax,xmeanmax,"log <f>","log #(<f>)",IOSTREAMV);
  IOSTREAMC <<IOSTREAMENDL<<" -> log-log histogram of the right tail of <f>:"<<IOSTREAMENDL;
  PlotLogLogHistogramInASCII(xmean,samples,x_ymax,xmeanmax,"log <f>","log #(<f>)",IOSTREAMV);
  
  // re-sample right tail of histogram (new method)
  double *x = NULL;
  double *y = NULL;
  double *w = NULL;
  double c0,c1,cov00,cov01,cov11,chisq;
  const long nr_observations = long(0.2*round(sqrt(samples)));
  Util_CreateFakeLogLogHistogram(&x,&y,&w,xmean,samples,x_ymax,xmeanmax,nr_observations);
  
  // make linear weighted fit using GSL
  gsl_fit_wlinear(x,1,w,1,y,1,nr_observations,&c0,&c1,&cov00,&cov01,&cov11,&chisq);
  IOSTREAMC <<" -> best fit: y(x) = "<<c0<<" + ("<<c1<<")*x; (chisq = "<<chisq<<")"<<IOSTREAMENDL;

  // IOSTREAMC <<"vectorxy = ";
  // Util_CoordinatedForMathematica(x,y,nr_observations,IOSTREAMV);
  // IOSTREAMC <<"vectorxw = ";
  // Util_CoordinatedForMathematica(x,w,nr_observations,IOSTREAMV);
  // IOSTREAMC <<"c0 = "<<c0<<";"<<IOSTREAMENDL;
  // IOSTREAMC <<"c1 = "<<c1<<";"<<IOSTREAMENDL<<IOSTREAMENDL;
  
  // tranforming everything back to linear space (where the line in log-log corresponds
  // to the power-law fit we wanted)
  IOSTREAMC <<" -> transforming back to linear space:"<<IOSTREAMENDL;
  for(int i=0; i<nr_observations; i++)
  {
    x[i] = exp(x[i]);
    y[i] = exp(y[i]);
    // w[i] = exp(w[i]); Unsinn, das ist ja normiert und so.

    // set up new weights in linear space
    w[i] = exp(-pow((x[i]-(x_ymax+xmeanmax)/2.)/((xmeanmax-x_ymax)/4.),2.));
  }
  // normalize weights (actually, this would not be necessary but is nicer maybe.)
  w[0] = 0.; // evil hack for the peak... :-(
  for(int i=0; i<nr_observations; i++)
    w[i] = w[i]/total(w,nr_observations);
    
  double coeffA = exp(c0);
  double coeffGamma = c1;
  IOSTREAMC <<" -> best fit in linear space: p(f) = "<<coeffA<<" * f^("<<coeffGamma<<")"<<IOSTREAMENDL;

  double xthresh;
  double* deviations = new double[nr_observations];
  for(int i=0; i<nr_observations; i++)
    deviations[i] = coeffA*pow(x[i],coeffGamma)-y[i];
  xthresh = 0.5*sqrt(gsl_stats_wvariance(w,1,deviations,1,nr_observations));
  IOSTREAMC <<" -> threshold for deviation: "<<xthresh<<IOSTREAMENDL;

  // IOSTREAMC <<" -> x = ";
  // display_subset(x,nr_observations,IOSTREAMV);
  // IOSTREAMC <<" -> w = ";
  // display_subset(w,nr_observations,IOSTREAMV);
  // IOSTREAMC <<" -> deviations = ";
  // display_subset(deviations,nr_observations,IOSTREAMV);

  // IOSTREAMC <<"debug: vectorxy = ";
  // Util_CoordinatedForMathematica(x,y,nr_observations,IOSTREAMV);
  // IOSTREAMC <<"debug: vectorxw = ";
  // Util_CoordinatedForMathematica(x,w,nr_observations,IOSTREAMV);
  IOSTREAMC <<"A = "<<coeffA<<";"<<IOSTREAMENDL;
  IOSTREAMC <<"gamma = "<<coeffGamma<<";"<<IOSTREAMENDL<<IOSTREAMENDL;

  // Von der Hälfte der Werte (auf lin. Skala) ab gehen wir nach unten (f=0), bis die
  // Abweichung vom Fit überschwellig wird:
  bool FoundConditioningLevel = false;
  int debug_count = 0;
  if(coeffGamma<0.)
  {
    for(int i=0; i<nr_observations; i++)
    {
      // find greatest x on left side for which the deviations are above theshold
      if(((x[i]-x_ymax)/(xmeanmax-x_ymax)<0.5) && (abs(deviations[i])>xthresh))
      {
        debug_count++;
        if(!FoundConditioningLevel)
        {
          FoundConditioningLevel = true;
          xresultlevel = x[i];
        }
        else if(x[i]>xresultlevel)
          xresultlevel = x[i];
      }
      else break;
    }
  }
  else IOSTREAMC <<" -> Warning: power law exponent is not negative as expected!"<<IOSTREAMENDL;
  IOSTREAMC <<" -> Debug: found "<<debug_count<<" points with superthreshold deviations."<<IOSTREAMENDL;
  
  if(!FoundConditioningLevel)
  {
    IOSTREAMC <<" -> Warning: Conditioning level could not be found, override at 1/3 between max and end!";
    IOSTREAMC <<IOSTREAMENDL;
    xresultlevel = 0.33333*(xmeanmax-x_ymax) + x_ymax;
  }
  
  // IOSTREAMC <<" -> number of observations for variance estimate: "<<nr_observations<<IOSTREAMENDL;
  // IOSTREAMC <<" -> std. deviation of difference to power-law fit (pseudo-weighted): ";
  // IOSTREAMC <<sqrt(gsl_stats_variance(deviations,1,nr_observations))<<IOSTREAMENDL;
  
  // free memory
  Util_FreeFakeHistogram(x,y,w);
  free_time_series_memory(xmean);
  delete[] deviations;
  return xresultlevel;
};

void Test_SetMinimaToZero(double** data, unsigned int size, long samples)
{
  double minimum;
  for(unsigned int i=0; i<size; i++)
  {
    minimum = smallest(data[i],samples);
    for(long t=0; t<samples; t++) data[i][t] -= minimum;
    assert(!(smallest(data[i],samples)<0.));
  }
};

void PlotHistogramInASCII(double* data, int samples, double range_min, double range_max, std::string xlabel, std::string ylabel, IOSTREAMH)
{
  PlotHistogramInASCII(false,false,data,samples,range_min,range_max,xlabel,ylabel,IOSTREAMV);
};
void PlotLogHistogramInASCII(double* data, int samples, double range_min, double range_max, std::string xlabel, std::string ylabel, IOSTREAMH)
{
  PlotHistogramInASCII(false,true,data,samples,range_min,range_max,xlabel,ylabel,IOSTREAMV);
};
void PlotLogLogHistogramInASCII(double* data, int samples, double range_min, double range_max, std::string xlabel, std::string ylabel, IOSTREAMH)
{
  PlotHistogramInASCII(true,true,data,samples,range_min,range_max,xlabel,ylabel,IOSTREAMV);
};
void PlotHistogramInASCII(bool xlogscaleQ, bool ylogscaleQ, double* data, int samples, double range_min, double range_max, std::string xlabel, std::string ylabel, IOSTREAMH)
{
  assert(range_min!=range_max);
  const int histo_bins = std::max(20,std::min(65,int(round(sqrt(samples)))));
  // IOSTREAMC <<" -> debug: histo_bins = "<<histo_bins<<IOSTREAMENDL;
  
  // create and fill histogram
  long* histo = new long[histo_bins];
  memset(histo, 0, histo_bins*sizeof(long));
  if(!xlogscaleQ)
    for(long t=0; t<samples; t++)
      histo[int(discretize(data[t],range_min,range_max,histo_bins))] += 1;
  else
    for(long t=0; t<samples; t++)
      histo[int(discretize(log(data[t]),log(range_min),log(range_max),histo_bins))] += 1;

  // convert histogram to floating point array (and free integer histogram)
  double* histoD = new double[histo_bins];
  for (int i=0; i<histo_bins; i++)
    histoD[i] = double(histo[i]);
  delete[] histo;

  if(ylogscaleQ)
    for(int i=0; i<histo_bins; i++)
    {
      if(histoD[i]>exp(1.)) histoD[i] = log(histoD[i]);
      else histoD[i] = 0.;
    }
  
  // find maximum of histogram
  double max_histo = 0.;
  for (unsigned int i=0; i<histo_bins; i++)
    if(histoD[i]>max_histo) max_histo = histoD[i];
  double min_histo = max_histo;
  for (unsigned int i=0; i<histo_bins; i++)
    if(histoD[i]<min_histo) min_histo = histoD[i];
       
  // draw histogram
  IOSTREAMC <<"^";
  if(ylogscaleQ) IOSTREAMC <<" (log)";
  IOSTREAMC <<IOSTREAMENDL;
  for(int line=0; line<HEIGHT_OF_ASCII_PLOTS; line++)
  {
    IOSTREAMC <<"|";
    for(int row=0; row<histo_bins; row++)
    {
      // if(histoD[row]/double(max_histo)>=(1.-double(line)/double(HEIGHT_OF_ASCII_PLOTS))) IOSTREAMC <<"#";
      if(histoD[row]/max_histo>=(1.-(double(line)+0.5)/double(HEIGHT_OF_ASCII_PLOTS)))
      {
        if(histoD[row]/max_histo>=(1.-(double(line))/double(HEIGHT_OF_ASCII_PLOTS)))
          IOSTREAMC <<":";
        else IOSTREAMC <<".";
      }
      else IOSTREAMC <<" ";
    }
    IOSTREAMC <<IOSTREAMENDL;
  }
  IOSTREAMC <<"+";
  for(int row=0; row<histo_bins; row++) IOSTREAMC <<"-";
  IOSTREAMC <<">";
  if(xlogscaleQ) IOSTREAMC <<" (log)";
  IOSTREAMC <<IOSTREAMENDL;
  
  // label axes
  IOSTREAMC <<"x-axis: "<<xlabel<<", ";
  IOSTREAMC <<"range "<<range_min<<" – "<<range_max<<IOSTREAMENDL;
  IOSTREAMC <<"y-axis: "<<ylabel<<", ";
  if(!ylogscaleQ) IOSTREAMC <<"range "<<min_histo<<" – "<<max_histo<<IOSTREAMENDL;
  else IOSTREAMC <<"range "<<long(exp(min_histo))<<" – "<<long(exp(max_histo))<<IOSTREAMENDL;
  
  // free memory
  delete[] histoD;
};

/* void PlotInASCII(double* data, int x_length, double range_min, double range_max, IOSTREAMD);
void PlotInASCII(double** data, int x_length, double range_min, double range_max, IOSTREAMD)
{
  assert(range_min!=range_max);
  //   const int histo_bins = std::max(20,std::min(65,int(round(sqrt(samples)))));
  //   // IOSTREAMC <<" -> debug: histo_bins = "<<histo_bins<<IOSTREAMENDL;
  //   
  //   // create and fill histogram
  //   long* histo = new long[histo_bins];
  //   memset(histo, 0, histo_bins*sizeof(long));
  //   if(!xlogscaleQ)
  //     for(long t=0; t<samples; t++)
  //       histo[int(discretize(data[t],range_min,range_max,histo_bins))] += 1;
  //   else
  //     for(long t=0; t<samples; t++)
  //       histo[int(discretize(log(data[t]),log(range_min),log(range_max),histo_bins))] += 1;
  // 
  //   // convert histogram to floating point array (and free integer histogram)
  //   double* histoD = new double[histo_bins];
  //   for (int i=0; i<histo_bins; i++)
  //     histoD[i] = double(histo[i]);
  //   delete[] histo;
  // 
  //   if(ylogscaleQ)
  //     for(int i=0; i<histo_bins; i++)
  //     {
  //       if(histoD[i]>exp(1.)) histoD[i] = log(histoD[i]);
  //       else histoD[i] = 0.;
  //     }
  //   
  //   // find maximum of histogram
  //   double max_histo = 0.;
  // for (unsigned int i=0; i<histo_bins; i++)
  //     if(histoD[i]>max_histo) max_histo = histoD[i];
  //   double min_histo = max_histo;
  // for (unsigned int i=0; i<histo_bins; i++)
  //     if(histoD[i]<min_histo) min_histo = histoD[i];
       
  // draw histogram
  IOSTREAMC <<"^";
  // if(ylogscaleQ) IOSTREAMC <<" (log)";
  IOSTREAMC <<IOSTREAMENDL;
  for(int line=0; line<HEIGHT_OF_ASCII_PLOTS; line++)
  {
    IOSTREAMC <<"|";
    for(int row=0; row<x_length; row++)
    {
      // if(histoD[row]/double(max_histo)>=(1.-double(line)/double(HEIGHT_OF_ASCII_PLOTS))) IOSTREAMC <<"#";
      if(histoD[row]/max_histo>=(1.-(double(line)+0.5)/double(HEIGHT_OF_ASCII_PLOTS)))
      {
        if(histoD[row]/max_histo>=(1.-(double(line))/double(HEIGHT_OF_ASCII_PLOTS)))
          IOSTREAMC <<":";
        else IOSTREAMC <<".";
      }
      else IOSTREAMC <<" ";
    }
    IOSTREAMC <<IOSTREAMENDL;
  }
  IOSTREAMC <<"+";
  for(int row=0; row<histo_bins; row++) IOSTREAMC <<"-";
  IOSTREAMC <<">";
  if(xlogscaleQ) IOSTREAMC <<" (log)";
  IOSTREAMC <<IOSTREAMENDL;
  
  // label axes
  IOSTREAMC <<"x-axis: "<<xlabel<<", ";
  IOSTREAMC <<"range "<<range_min<<" – "<<range_max<<IOSTREAMENDL;
  IOSTREAMC <<"y-axis: "<<ylabel<<", ";
  if(!ylogscaleQ) IOSTREAMC <<"range "<<min_histo<<" – "<<max_histo<<IOSTREAMENDL;
  else IOSTREAMC <<"range "<<long(exp(min_histo))<<" – "<<long(exp(max_histo))<<IOSTREAMENDL;
  
  // free memory
  delete[] histoD;
}; */



double Util_FindPeakInHistogram(const double* data, const long samples, const double range_min, const double range_max, const int histo_bins)
{
  // create and fill histogram
  long* histo = new long[histo_bins];
  memset(histo, 0, histo_bins*sizeof(long));
  for(long t=0; t<samples; t++)
    histo[int(discretize(data[t],range_min,range_max,histo_bins))] += 1;

  double x = range_min;
  long y = 0;
  
  for(int i=0; i<histo_bins; i++)
  {
    if(histo[i]>y)
    {
      // found a higher point:
      x = (double(i)+0.5)/double(histo_bins)*(range_max-range_min) + range_min;
      y = histo[i];
    }
  }
  
  delete[] histo;
  return x;
};
void Util_CreateFakeLogLogHistogram(double** x, double** y, double** w, const double* data, const long samples, const double range_min, const double range_max, const int histo_bins)
{
  // create and fill histogram
  long* histo = new long[histo_bins];
  memset(histo, 0, histo_bins*sizeof(long));
  for(long t=0; t<samples; t++)
    histo[int(discretize(log(data[t]),log(range_min),log(range_max),histo_bins))] += 1;

  *x = new double[histo_bins];
  *y = new double[histo_bins];
  *w = new double[histo_bins];
  for(int i=0; i<histo_bins; i++)
  {
    // set up x
    if(histo[i]>2) (*y)[i] = log(double(histo[i]));
    else (*y)[i] = 0.;
    
    // set up x and w
    (*x)[i] = log(double(i+1)/double(histo_bins)*(range_max-range_min) + range_min);
    (*w)[i] = exp(-pow((double(i)-double(histo_bins)/2.)/(double(histo_bins)/6.),2.));
  }
  // normalize weights
  for(int i=0; i<histo_bins; i++)
    (*w)[i] = ((*w)[i])/total(*w,histo_bins);
  
  delete[] histo;
};
void Util_FreeFakeHistogram(double* x, double* y, double* w)
{
  delete[] x;
  delete[] y;
  delete[] w;
};

void Util_CoordinatedForMathematica(double*x, double*y, int length, IOSTREAMH)
{
  IOSTREAMC <<"{";
  for(int i=0; i<length; i++)
  {
    if(i>0) IOSTREAMC <<",";
    IOSTREAMC <<"{"<<x[i]<<","<<y[i]<<"}";
  }
  IOSTREAMC <<"};"<<IOSTREAMENDL;
};

double** generate_conditioned_time_series_by_glueing(double** data, const int size, double* xmean, const long StartSampleIndex, const long EndSampleIndex, const double condlevel, unsigned long* available_samples, IOSTREAMH)
{
  // determine number of samples that will be available
  *available_samples = 0;
  for(long t=StartSampleIndex; t<=EndSampleIndex; t++)
  {
    if(xmean[t]<condlevel) *available_samples += 1;
  }
  
  // allocate memory
  double** result = new double*[size];
  for(int i=0; i<size; i++)
    result[i] = new double[*available_samples];
  
  // generate glued signal
  unsigned long added_samples = 0;
  for(long t=StartSampleIndex; t<=EndSampleIndex; t++)
  {
    if(xmean[t]<condlevel)
    {
      for(int i=0; i<size; i++)
        result[i][added_samples] = data[i][t];
      added_samples++;
    }
  }
  assert(added_samples==*available_samples);
  
  
  IOSTREAMC <<" -> global conditioning level "<<condlevel<<": ";
  IOSTREAMC <<((100.*added_samples)/(EndSampleIndex-StartSampleIndex));
  IOSTREAMC <<"% are below threshold. "<<IOSTREAMENDL;
  
  return result;
};

// code taken from:
// http://code.google.com/p/yaml-cpp/wiki/HowToParseADocument
double** read_positions_from_YAML(std::string YAMLfilename, unsigned int size, IOSTREAMH)
{
#ifndef ENABLE_YAML_IMPORT_AT_COMPILE_TIME
  IOSTREAMC <<"error: YAML input disabled at compile time!"<<IOSTREAMENDL;
#else
  
  YAML::Node yamldoc;
  std::ifstream fin(YAMLfilename.c_str());
  if(!fin.is_open())
  {
    IOSTREAMC <<"error: YAML input file '"<<YAMLfilename<<"' not found!"<<IOSTREAMENDL;
    exit(1);
  }
  else IOSTREAMC <<"-> loading YAML input file '"<<YAMLfilename<<"' ..."<<IOSTREAMENDL;

  double** positions = NULL;
  unsigned int nr_of_position_entries = 0;
  try {
    YAML::Parser parser(fin);
    parser.GetNextDocument(yamldoc);
    
    // test if 'size' tag is present and if entry matches control file
    std::string name;
    yamldoc["size"] >> name;
    unsigned int size_read = atoi(name.c_str());
    // IOSTREAMC <<"loading from YAML file: size = "<<size_read<<IOSTREAMENDL;
    if(size!=size_read) {
      IOSTREAMC <<"error while loading from YAML file: key 'size' does not match size parameter.";
      IOSTREAMC <<IOSTREAMENDL;
      exit(1);
    }
    
    // iterate through nodes and extract positions
    positions = new double*[size];
    for(unsigned int i=0; i<size; i++)
      positions[i] = new double[2];
    if(const YAML::Node *Nodes = yamldoc.FindValue("nodes")) {
      // make sure we are at the right place
      // IOSTREAMC <<"debug: found record of "<<(*Nodes).size()<<" nodes."<<IOSTREAMENDL;
      // assert((**Nodes).getType()==YAML::CT_SEQUENCE);
      assert((*Nodes).size()==size);
      for(unsigned int i=0; i<size; i++)
      {
        const YAML::Node& myNode = (*Nodes)[i];

        // 1.) read id
        std::string id;
        *myNode.FindValue("id") >> id;
        unsigned int read_id = atoi(id.c_str());
        // IOSTREAMC <<"debug: node #"<<read_id<<IOSTREAMENDL;

        // 2.) read position
        if(const YAML::Node *myPos = myNode.FindValue("pos")) {
          nr_of_position_entries++;
          // IOSTREAMC <<"debug: found position entry for node #"<<read_id<<"."<<IOSTREAMENDL;
          std::string readfloat;
          (*myPos)[0] >> readfloat;
          positions[read_id-1][0] = atof(readfloat.c_str());
          (*myPos)[1] >> readfloat;
          positions[read_id-1][1] = atof(readfloat.c_str());

          // IOSTREAMC <<"debug: found position of node #"<<read_id<<": ";
          // IOSTREAMC <<positions[read_id-1][0]<<", "<<positions[read_id-1][1]<<IOSTREAMENDL;
        }
        else {
          IOSTREAMC <<"error while loading from YAML file: could not find position entry for node #";
          IOSTREAMC <<read_id<<"."<<IOSTREAMENDL;
          exit(1);
        }
      }
    }  else {
      IOSTREAMC <<"error while loading from YAML file: key 'nodes' does not exist."<<IOSTREAMENDL;
      exit(1);
    };

  }
  catch(YAML::ParserException& e) {
    IOSTREAMC << e.what() << "\n";
    exit(1);
  }
  
  IOSTREAMC <<"-> position entries for "<<nr_of_position_entries<<" nodes have been loaded."<<IOSTREAMENDL;
  fin.close();
  return positions;
#endif
};
void free_position_memory(double** pos, unsigned int size) {
  free_time_series_memory(pos,size);
};

double norm(double* pointA, double* pointB) {
  return(sqrt(pow(pointA[0]-pointB[0],2.)+pow(pointA[1]-pointB[1],2.)));
};
double** clone_time_series(double** data, unsigned int size, long samples) {
  double **data_copy = new double*[size];
  for(unsigned int i=0; i<size; i++)
  {
    data_copy[i] = new double[samples];
    memcpy(data_copy[i],data[i],samples*sizeof(double));
  }
  return data_copy;
};
void apply_light_scattering_to_time_series(double** data, unsigned int size, long samples, std::string YAMLfilename, double sigma_scatter, double amplitude_scatter, IOSTREAMH)
{
  // clone data memory
  double **data_copy = clone_time_series(data,size,samples);
  
  // read node positions from YAML
  double** positions = read_positions_from_YAML(YAMLfilename,size,IOSTREAMV);
  double dist;
  double* ScatterAmplitudes = new double[size];
  
  // std::cout <<"-> running "<<std::flush; 
  for(unsigned int i=0; i<size; i++)
  {
    // std::cout <<"."<<std::flush;
     // this assumes that light scattering can be filtered out to some extent (other nodes are
     // multiplied by scatter_amplitude)
    ScatterAmplitudes[i] = 0.;
    for(unsigned int j=0; j<size; j++)
      // calculate the amount of scattering that (j) has on (i)
      if(j!=i) {
        dist = norm(positions[i],positions[j]);
        ScatterAmplitudes[j] = amplitude_scatter*exp(-pow(dist/sigma_scatter,2.));
        assert(ScatterAmplitudes[j]>=0.);
      }

    // apply scattering
    for(long t=0; t<samples; t++)
      for(unsigned int j=0; j<size; j++)
        data[i][t] +=  ScatterAmplitudes[j]*data_copy[j][t]; // this would be buggy if the amplitudes could be negative...
  }
  // std::cout <<std::endl;
  // de-allocate memory
  delete[] ScatterAmplitudes;
  free_time_series_memory(data_copy,size);
  free_position_memory(positions,size);
};

double SphericalUnitSurface(int r) {
  return (r*pow(PI,r/2.))/(gsl_sf_gamma(double(r)/2.+1));
};
double gsl_norm(const gsl_vector* vecA, const gsl_vector* vecB, int dim) {
  double result = 0.;
  for(int i=0; i<dim; i++)
    result += pow(gsl_vector_get(vecA,i)-gsl_vector_get(vecB,i),2.);
  return sqrt(result);
};
double gsl_quicknorm(const gsl_vector* vecA, const gsl_vector* vecB, int dim, double bound) {
  double result = 0.;
  double OneDdist;
  for(int i=0; i<dim; i++) {
    OneDdist = abs(gsl_vector_get(vecA,i)-gsl_vector_get(vecB,i));
    if (OneDdist > bound) return 2.*bound; // break if distance in 1D is already to large
    result += OneDdist*OneDdist;
  }
  return sqrt(result);
};

long double DifferentialEntropy(gsl_vector** data, const int dim, const long samples)
{
  // reference:
  // Victor. Binless strategies for estimation of information from neural data. Physical
  // Review E (2002) vol. 66 (5) pp. 51903: see there eq. 10
  long double Hdiff = 0.;
  double lowest_distance, distance_here;
  
  // cout <<"direkt: dim = "<<dim<<", samples = "<<samples<<endl;
  // find nearest neighbor distances (1st term in Hdiff)
  for(long s=0; s<samples; s++)
  {
    lowest_distance = std::numeric_limits<double>::max();
    // find sample closest to sample with index s
    for(long j=0; j<samples; j++) {
      if (s!=j) {
        distance_here = gsl_quicknorm(data[s],data[j],dim,lowest_distance);
        if(distance_here<lowest_distance) {
          lowest_distance = distance_here;
        }
      }
    }
    // lowest_distance = sqrt(lowest_distance);
    // std::cout <<"debug: s="<<s<<", lowest_distance="<<lowest_distance;
    // std::cout <<", distance_here="<<distance_here<<std::endl;

    Hdiff += log(lowest_distance);
    
    // cout <<"direkt: NN distance of sample #"<<s<<": "<<lowest_distance<<endl;
  }

  // std::cout <<"debug: Hdiff_sumonly = "<<Hdiff<<std::endl;
  Hdiff *= double(dim)/(double(samples)*log(2.));
  // std::cout <<"debug: Hdiff_1 = "<<Hdiff<<std::endl;
  
  // second term
  Hdiff += double(dim)*log(SphericalUnitSurface(dim)*double(samples-1)/double(dim))/log(2.);
  // std::cout <<"debug: Hdiff_12 = "<<Hdiff<<std::endl;

  // third term
  Hdiff += double(dim)*EULERGAMMA/log(2.);
  // std::cout <<"debug: Hdiff_123 = "<<Hdiff<<std::endl;
  
  return Hdiff;
};

long* generate_random_permutation(long samples, gsl_rng* GSLrandom)
{
  // initialize random number generator
  // gsl_rng_env_setup();
  // gsl_rng* GSLrandom = gsl_rng_alloc(GSL_RANDOM_NUMBER_GENERATOR);
  // gsl_rng_set(GSLrandom, 1234);

  long* shuffle_permutation = new long[samples];
  // generate random permutation (once for all)
  for(long t=0; t<samples; t++) shuffle_permutation[t] = t;
  gsl_ran_shuffle(GSLrandom,shuffle_permutation,samples,sizeof(long));
  
  // gsl_rng_free(GSLrandom);
  return shuffle_permutation;
};

long* generate_random_permutation(long samples, rawdata globalbins, unsigned long* AvailableSamples, long StartSampleIndex, long EndSampleIndex, rawdata* xglobal, gsl_rng* GSLrandom)
{
  if (globalbins<2) return generate_random_permutation(samples,GSLrandom);
  
  // initialize random number generator
  // gsl_rng_env_setup();
  // gsl_rng* GSLrandom = gsl_rng_alloc(GSL_RANDOM_NUMBER_GENERATOR);
  // gsl_rng_set(GSLrandom, 1234);
  
  long* shuffle_permutation = new long[samples];
  long t_condindex;

  for(rawdata g=0; g<globalbins; g++) {
    // cout <<"debug: g = "<<int(g)<<endl;
    if(AvailableSamples[g]>0) {
      long* conditioned_shuffle = new long[AvailableSamples[g]];
      // 1.) prepare shuffles sample indices coming from one particular globalbin
      t_condindex = 0;
      for(long t=StartSampleIndex; t<=EndSampleIndex; t++) {
        if(xglobal[t]==g) {
          conditioned_shuffle[t_condindex] = t;
          t_condindex++;
        }
      }
      // cout <<"t_condindex = "<<t_condindex<<", AvailableSamples[g] = "<<AvailableSamples[g]<<endl;
      assert(t_condindex==AvailableSamples[g]);
      // 2.) shuffle all those indices
      gsl_ran_shuffle(GSLrandom,conditioned_shuffle,AvailableSamples[g],sizeof(long));
      // 3.) insert indices at places where the globalbin has the correct value
      t_condindex = 0;
      for(long t=StartSampleIndex; t<=EndSampleIndex; t++) {
        if(xglobal[t]==g) {
          shuffle_permutation[t] = conditioned_shuffle[t_condindex];
          t_condindex++;
        }
      }
      assert(t_condindex==AvailableSamples[g]);
      delete[] conditioned_shuffle;
    }
  }

  // gsl_rng_free(GSLrandom);
  return shuffle_permutation;
};

void random_permutation(long** data, const long samples, gsl_rng* GSLrandom)
{
  // initialize random number generator
  // gsl_rng_env_setup();
  // gsl_rng* GSLrandom = gsl_rng_alloc(GSL_RANDOM_NUMBER_GENERATOR);
  // gsl_rng_set(GSLrandom, rand());

  // test if input data is valid (memchecker)
  for(long t=0; t<samples; t++) assert((*data)[t]>=0);
  gsl_ran_shuffle(GSLrandom,*data,samples,sizeof(long));
  
  // gsl_rng_free(GSLrandom);  
};

void geometric_permutation(long** data, const long samples, const long AutoCorrLength, gsl_rng* GSLrandom)
{
  for(long t=0; t<samples; t++) assert((*data)[t]>=0);

  const long LEGObits = samples/AutoCorrLength; // implicit round down
  // std::cout <<"debug: samples = "<<samples<<", AutoCorrLength = "<<AutoCorrLength<<", LEGObits = "<<LEGObits<<std::endl;
  if(LEGObits==0) return;

  long* LEGObits_permutation = new long[LEGObits+1];
  for(long t=0; t<LEGObits+1; t++) LEGObits_permutation[t] = t;
  // shuffle LEGObits (except for the last one)
  random_permutation(&LEGObits_permutation,LEGObits+0,GSLrandom);
  
  // commit
  long* data_copy = new long[samples];
  // memcpy(data_copy,*data,samples*sizeof(long)); short-hand
  // explicit:
  for(long t=0; t<samples; t++) {
    data_copy[t] = (*data)[t];
  }
  
  // for asserting that every value has been set
  for(long t=0; t<samples; t++) (*data)[t] = -1;
  long bitindex, loopindex;
  // for(long t=0; t<samples; t++) {
  //   bitindex = long(floor(double(t)/double(AutoCorrLength)));
  //   // bitindex = t/AutoCorrLength;
  //   assert(bitindex<=LEGObits);
  //   loopindex = t % AutoCorrLength;
  //   (*data)[t] = data_copy[loopindex + LEGObits_permutation[bitindex]*AutoCorrLength];
  //   // std::cout <<"geometric_permutation: t="<<t<<": permutation[t] = "<<data[t]<<std::endl;
  // }
  long t, t_copy;
  for(long bit=0; bit<LEGObits; bit++) {
    for(long bit_part=0; bit_part<AutoCorrLength; bit_part++) {
      t = bit_part + bit*AutoCorrLength;
      t_copy = bit_part + LEGObits_permutation[bit]*AutoCorrLength;
      // std::cout <<"debug: t = "<<t<<", t_copy = "<<t_copy<<" (bit = "<<bit<<", bit_part = "<<bit_part<<")"<<std::endl;
      (*data)[t] = data_copy[t_copy];
    }
  }
  // std::cout <<"--- fill in rest ---"<<std::endl;
  for(long t_rest = t+1; t_rest<samples; t_rest++) {
    // std::cout <<"t_rest = "<<t_rest<<", LEGObits="<<LEGObits<<", samples = "<<samples<<std::endl;
    (*data)[t_rest] = data_copy[t_rest];
  }
  for(long t=0; t<samples; t++) assert((*data)[t]>=0);
  
  delete[] data_copy;
  delete[] LEGObits_permutation;
};

long* generate_random_geometric_permutation(long samples, rawdata globalbins, rawdata* xglobal, long AutoCorrLength, gsl_rng* GSLrandom)
{
  long* shuffle_permutation = new long[samples];
  for(long t=0; t<samples; t++) shuffle_permutation[t] = -1; // to test if all values have been set
  long* shuffle_permutation_per_globalbin;
  long t_sample;
  
  // new counter (ignoring Start and EndSampleIndex)
  unsigned long* AllAvailableSamples = new unsigned long[globalbins];
  for(rawdata g=0; g<globalbins; g++) AllAvailableSamples[g] = 0;
  for(long t=0; t<samples; t++) AllAvailableSamples[long(xglobal[t])] += 1;
  
  for(rawdata g=0; g<globalbins; g++)
  {
    // std::cout <<"generate_random_geometric_permutation: starting g="<<int(g)<<" ..."<<std::endl;
    // LEGObits = long(floor(double(AllAvailableSamples[g])/double(AutoCorrLength)));
    // LEGObits_permutation = new long[LEGObits+1];
    
    // 1.) List mit Samples in gbin erstellen
    shuffle_permutation_per_globalbin = new long[AllAvailableSamples[g]];
    t_sample = 0;
    for(long t=0; t<samples; t++) {
      if(xglobal[t]==g) {
        shuffle_permutation_per_globalbin[t_sample] = t;
        t_sample++;
        // std::cout <<"g="<<int(g)<<": (t="<<t<<", t_sample="<<t_sample<<")\t"<<std::flush;
      }
    }
    assert(t_sample == AllAvailableSamples[g]);
    // std::cout <<std::endl<<"-> set up shuffle_permutation_per_globalbin."<<std::endl;
        
    // 2.) Die Liste von (1.) geometrisch Shuffeln
    geometric_permutation(&shuffle_permutation_per_globalbin,AllAvailableSamples[g],AutoCorrLength,GSLrandom);
    
    // 3.) Comitten pro Globalbin
    long t_in_my_gbin = 0;
    for(long t=0; t<samples; t++) {
      if(xglobal[t]==g) {
        shuffle_permutation[t] = shuffle_permutation_per_globalbin[t_in_my_gbin];
        t_in_my_gbin++;
      }
    }
    assert(t_in_my_gbin == AllAvailableSamples[g]);
    
    delete[] shuffle_permutation_per_globalbin;
    // std::cout <<"generate_random_geometric_permutation: g="<<int(g)<<" done."<<std::endl;
  }
  
  return shuffle_permutation;
};

double AutoCorrelation(double* data, const long samples, const long lag, bool Abs)
{
  const long effective_samples = samples-2*abs(lag);

  double* shifted_copy = new double[effective_samples];
  for(long s=0; s<effective_samples; s++) {
    shifted_copy[s] = data[s+lag];
  }
  double acc = gsl_stats_correlation(data,1,shifted_copy,1,effective_samples);
  if(Abs && acc<0.) acc *= -1;
  
  delete[] shifted_copy;
  // std::cout <<"debug: AutoCorrelation: acc = "<<acc<<std::endl;
  return acc;
};
double AutoCorrelationTimeScale(double* data, const long samples, const long max_lag, IOSTREAMH)
{
  double* x_coords = new double[max_lag+1];
  x_coords[0] = 0.;
  double* y_coords = new double[max_lag+1];
  y_coords[0] = log(1.);
  // ...because the ACC for lag=0 is trivially 1.
  
  for(long lag=1; lag<=max_lag; lag++) {
    x_coords[lag] = lag;
    y_coords[lag] = log(AutoCorrelation(data,samples,lag,true));
    // std::cout <<"debug: AutoCorrelationTimeScale: added coordinate ("<<x_coords[lag]<<", "<<y_coords[lag]<<")"<<std::endl;
  }
  // std::cout <<"debug: AutoCorrelationTimeScale: in Mathematica format:"<<std::endl;
  Util_CoordinatedForMathematica(x_coords,y_coords,max_lag+1,IOSTREAMV);
  
  double c0,c1,cov00,cov01,cov11,sumsq;
  gsl_fit_linear(x_coords,1,y_coords,1,max_lag+1,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
  
  return -1./c1;
};


void write_result(double **array, long size, std::string outputfile_results_name, IOSTREAMH)
{
  char* name = new char[outputfile_results_name.length()+1];
  strcpy(name,outputfile_results_name.c_str());
  ofstream fileout1(name);
  delete[] name;
  if (fileout1 == NULL) {
    IOSTREAMC <<IOSTREAMENDL<<"error in write_result: cannot open output file!"<<IOSTREAMENDL;
    exit(1);
  }   

  fileout1.precision(OUTPUTNUMBER_PRECISION);
  fileout1 <<fixed;
  fileout1 <<"{";
  for(unsigned int j=0; j<size; j++) {
    if(j>0) fileout1<<",";
    fileout1 <<"{";
    for(unsigned int i=0; i<size; i++) {
      if (i>0) fileout1<<",";
      fileout1 <<(double)array[j][i];
    }
    fileout1 <<"}"<<endl;
  }
  fileout1 <<"}"<<endl;

  IOSTREAMC <<"result saved to file."<<IOSTREAMENDL;
};


void write_multidim_result(double ***array, unsigned int dimens, long size, std::string outputfile_results_name, IOSTREAMH)
{
  char* name = new char[outputfile_results_name.length()+1];
  strcpy(name,outputfile_results_name.c_str());
  ofstream fileout1(name);
  delete[] name;
  if (fileout1 == NULL) {
    IOSTREAMC <<IOSTREAMENDL<<"error in write_multidim_result: cannot open output file!"<<IOSTREAMENDL;
    exit(1);
  }   

  fileout1.precision(OUTPUTNUMBER_PRECISION);
  fileout1 <<fixed;
  fileout1 <<"{";
  for(unsigned int j=0; j<size; j++) {
    if(j>0) fileout1<<",";
    fileout1 <<"{";
    for(unsigned int i=0; i<size; i++) {
      if (i>0) fileout1<<",";
      fileout1 <<"{";
      for(int k=0; k<dimens; k++) {
        if (k>0) fileout1<<",";
        fileout1 <<array[j][i][k];
      }
      fileout1 <<"}";
    }
    fileout1 <<"}"<<endl;
  }
  fileout1 <<"}"<<endl;

  IOSTREAMC <<"result saved to file."<<IOSTREAMENDL;
};


std::string bool2textMX(bool value)
{
  if (value) return "True";
  else return "False";
};

void apply_baseline_correction(double* data, long samples)
{
  long i_start, i_end, i_width;
  double* temp_time_series = new double[samples];
  // copy one time series to temp
  memcpy(temp_time_series,data,samples*sizeof(double));

#ifndef SPEEDUP_BASELINE_CORRECTION
  long double temp_mean;

  for(long tt=0; tt<samples; tt++) {
    i_start = std::max((long)0,tt-BASELINE_CORRECTION_BANDWIDTH);
    i_end = std::min(samples-1,tt+BASELINE_CORRECTION_BANDWIDTH);
    i_width = i_end-i_start+1;
    assert(i_width>0);
    
    // calulate mean signal in this region
    temp_mean = 0.;
    for(long tt2=i_start; tt2<=i_end; tt2++) {
      temp_mean += temp_time_series[tt2];
    }
    temp_mean /= (long double)i_width;
    
    // apply correction
    data[tt] -= (double)temp_mean;
    // cout <<"debug baseline corr.: xresult["<<i<<"]["<<tt<<"]: pre = "<<temp_time_series[tt]<<", post = "<<xresult[i][tt]<<endl;
  }
#else
  long double moving_window_total = 0.;
  long actual_start = 0;
  long actual_end = -1;

  for(long tt=0; tt<samples; tt++) {
    i_start = std::max((long)0,tt-BASELINE_CORRECTION_BANDWIDTH);
    i_end = std::min(samples-1,tt+BASELINE_CORRECTION_BANDWIDTH);
    i_width = i_end-i_start+1;
    assert(i_width>0);
    
    // shift around moving window
    while(actual_end < i_end) {
      moving_window_total += (long double)(temp_time_series[actual_end+1]);
      actual_end++;
    }
    while(actual_start < i_start) {
      moving_window_total -= (long double)(temp_time_series[actual_start]);
      actual_start++;
      // std::cout <<"debug: reduction of "<<(long double)(temp_time_series[actual_start])<<"!"<<std::endl;
    }
    assert(i_end == actual_end);
    assert(i_start == actual_start);
    
    // apply correction
    data[tt] -= (double)(moving_window_total/((long double)i_width));
    
    // std::cout <<"debug: tt="<<tt<<": i_start="<<i_start<<", i_end="<<i_end<<", moving_window_total="<<moving_window_total<<", data[tt]="<<data[tt]<<std::endl;
  }
#endif

  delete[] temp_time_series;
};
