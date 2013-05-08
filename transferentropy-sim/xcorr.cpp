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
// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
// #include <gsl/gsl_vector.h>

#include "../olav.h"
#include <sim_main.h>
#include "../te-datainit.h"

// #ifndef INLINE
// #define INLINE extern inline
// #endif

#define REPORTS 25
// #define SHOW_DETAILED_PROGRESS

#undef SEPARATED_OUTPUT

// #define GSL_RANDOM_NUMBER_GENERATOR gsl_rng_default
#define GSL_RANDOM_NUMBER_GENERATOR gsl_rng_ranlxs2

// #undef ENABLE_ADAPTIVE_BINNING_AT_COMPILE_TIME

#define CROSSCORRELATION_MAX_LAG 4

using namespace std;

// typedef unsigned char rawdata;

time_t start, middle, end, now;

class Kernel;

int main(int argc, char* argv[])
{
  SimControl<Kernel> simc;
  time(&start);
  int ret = simc.simulate(argc, argv);
  return 0;
};

class Kernel
{
public:
  unsigned long iteration;
  unsigned int size;
  // unsigned int bins, globalbins;
  unsigned int globalbins;
  // unsigned long mag der SimKernel irgendwie nicht?
  long samples;
  long StartSampleIndex, EndSampleIndex;
  bool EqualSampleNumberQ;
  long MaxSampleNumberPerBin;
  // unsigned long * AvailableSamples;
  unsigned long AvailableSamples;
  unsigned int word_length;
  double std_noise;
  string inputfile_name;
  string outputfile_results_name;
  string outputfile_pars_name;
  string spikeindexfile_name, spiketimesfile_name;
  string FluorescenceModel;

  double input_scaling;
  double cutoff;
  double tauF;
  double tauCa;
  double fluorescence_saturation;
  double DeltaCalciumOnAP;
  
  // parameters for light scattering
  std::string YAMLfilename;
  double SigmaScatter;
  double AmplitudeScatter;
  
  bool OverrideRescalingQ; // if, for example the input are pure spike data (integers)
  bool HighPassFilterQ; // actually, this transforms the signal into the difference signal
  bool InstantFeedbackTermQ;
#ifdef ENABLE_ADAPTIVE_BINNING_AT_COMPILE_TIME
  bool AdaptiveBinningQ; // rubbish
#endif
  bool IncludeGlobalSignalQ;
  bool GenerateGlobalFromFilteredDataQ;
  double GlobalConditioningLevel;
  bool RelativeGlobalConditioningLevelQ;
  
  bool ContinueOnErrorQ;
  bool skip_the_rest;
  
  // bool AutoBinNumberQ;
  bool AutoConditioningLevelQ;
  
  gsl_rng* GSLrandom;

#ifndef SEPARATED_OUTPUT
  double **xresult;
#else
  double ***xresult;
#endif
  double** xdatadouble;
  double* xglobaldouble;
  double** xdatadoubleglue;

  void initialize(Sim& sim)
  {
    iteration = sim.iteration();
    sim.io <<"Init: iteration "<<iteration<<", process "<< sim.process()<<Endl;
    time(&now);
    sim.io <<"time: ";
    sim.io <<"elapsed "<<sec2string(difftime(now,start));
    sim.io <<", ETA "<<ETAstring(sim.iteration()-1,sim.n_iterations(),difftime(now,start))<<Endl;
    
    // read parameters from control file
    sim.get("size",size);
    // bins = 0;
    // sim.get("AutoBinNumberQ",AutoBinNumberQ,false);
    // if(!AutoBinNumberQ) sim.get("bins",bins);
    sim.get("globalbins",globalbins,1);
    sim.get("samples",samples);
    sim.get("StartSampleIndex",StartSampleIndex,1);
    sim.get("EndSampleIndex",EndSampleIndex,samples-1);
    sim.get("EqualSampleNumberQ",EqualSampleNumberQ,false);
    sim.get("MaxSampleNumberPerBin",MaxSampleNumberPerBin,-1);

    sim.get("noise",std_noise,-1.);
    sim.get("appliedscaling",input_scaling,1.);
    sim.get("cutoff",cutoff,-1.);
    sim.get("tauF",tauF,-1.);
    sim.get("OverrideRescalingQ",OverrideRescalingQ,false);
    sim.get("HighPassFilterQ",HighPassFilterQ,false);
    sim.get("InstantFeedbackTermQ",InstantFeedbackTermQ,false);
    sim.get("saturation",fluorescence_saturation,-1.);
    sim.get("IncludeGlobalSignalQ",IncludeGlobalSignalQ,false);
    assert(IncludeGlobalSignalQ ^ (globalbins==1));
    sim.get("GenerateGlobalFromFilteredDataQ",GenerateGlobalFromFilteredDataQ,false);
    sim.get("AutoConditioningLevelQ",AutoConditioningLevelQ,false);
    if(!AutoConditioningLevelQ)
    {
      sim.get("GlobalConditioningLevel",GlobalConditioningLevel,-1.);
      if (GlobalConditioningLevel>0) assert(globalbins==2); // for now, to keep it simple
    }
    else GlobalConditioningLevel = -1.;
    sim.get("RelativeGlobalConditioningLevelQ",RelativeGlobalConditioningLevelQ,false);
    
    sim.get("inputfile",inputfile_name,"");
    sim.get("outputfile",outputfile_results_name);
    sim.get("outputparsfile",outputfile_pars_name);
    sim.get("spikeindexfile",spikeindexfile_name,"");
    sim.get("spiketimesfile",spiketimesfile_name,"");
    // make sure we either have a fluorescence input or both spike inputs
    if(!((inputfile_name!="") ^ ((spikeindexfile_name!="")&&(spiketimesfile_name!=""))))
    {
      sim.io <<"Error: Based on the parameters, it is not clear where your data should come from."<<Endl;
      exit(1);
    }
    sim.get("FluorescenceModel",FluorescenceModel,"");
    sim.get("DeltaCalciumOnAP",DeltaCalciumOnAP,50);
    sim.get("tauCa",tauCa,1000);

    // parameters for light scattering
    sim.get("YAMLfile",YAMLfilename,"");
    sim.get("SigmaScatter",SigmaScatter,-1.);
    sim.get("AmplitudeScatter",AmplitudeScatter,-1.);

    sim.get("ContinueOnErrorQ",ContinueOnErrorQ,false);

    // initialize random number generator
    gsl_rng_env_setup();
    GSLrandom = gsl_rng_alloc(GSL_RANDOM_NUMBER_GENERATOR);
    gsl_rng_set(GSLrandom, 1234);

    AvailableSamples = 0;
  };

  void execute(Sim& sim)
  {
    sim.io <<"------ xcorrelation-sim:v2 ------ olav, Wed 18 May 2011 ------"<<Endl;
    // time_t start, middle, end;

    sim.io <<"output file: "<<outputfile_results_name<<Endl;
    // Gespeichert wird später - hier nur Test, ob das Zielverzeichnis existiert
    write_parameters();
    skip_the_rest = false;

    sim.io <<"allocating memory..."<<Endl;
    try {
      xdatadouble = NULL;
      xglobaldouble = NULL;
      xdatadoubleglue = NULL;

      xresult = NULL;
#ifndef SEPARATED_OUTPUT
      xresult = new double*[size];
#else
      xresult = new double**[size];
#endif
      for(int i=0; i<size; i++)
      {
        // xdata[i] = new rawdata[samples];
        // memset(xdata[i], 0, samples*sizeof(rawdata));
#ifndef SEPARATED_OUTPUT
        xresult[i] = new double[size];
        memset(xresult[i], 0, size*sizeof(double));
#else
        xresult[i] = new double*[size];
        for(int i2=0; i2<size; i2++)
        {
          xresult[i][i2] = new double[globalbins];
          memset(xresult[i][i2], 0, globalbins*sizeof(double));
        }
#endif
      }
    
      // AvailableSamples = new unsigned long[2];
      sim.io <<" -> done."<<Endl;
            
      if(inputfile_name!="")
        sim.io <<"input file: \""<<inputfile_name<<"\""<<Endl;
      else {
        sim.io <<"input files:"<<Endl;
        sim.io <<"- spike indices: \""<<spikeindexfile_name<<"\""<<Endl;
        sim.io <<"- spike times: \""<<spiketimesfile_name<<"\""<<Endl;
      }
      
      if(inputfile_name=="") {
        sim.io <<"loading data and generating time series from spike data..."<<Endl;
        xdatadouble = generate_time_series_from_spike_data(spiketimesfile_name, spikeindexfile_name, size, int(round(tauF)), samples, FluorescenceModel, std_noise, fluorescence_saturation, cutoff, DeltaCalciumOnAP, tauCa, GSLrandom, sim);
      }
      else {
        sim.io <<"loading data from binary file..."<<Endl;
        xdatadouble = load_time_series_from_file(inputfile_name, size, samples, input_scaling, OverrideRescalingQ, std_noise, fluorescence_saturation, cutoff, GSLrandom, sim);
      }
      sim.io <<" -> done."<<Endl;

      if(AmplitudeScatter>0.) {
        sim.io <<"simulating light scattering..."<<Endl;
        apply_light_scattering_to_time_series(xdatadouble, size, samples, YAMLfilename, SigmaScatter, AmplitudeScatter, sim);
        sim.io <<" -> done."<<Endl;
      }
      
      sim.io <<"histogram of averaged signal:"<<Endl;
      double* xmean = generate_mean_time_series(xdatadouble,size,samples);
      PlotLogHistogramInASCII(xmean,samples,smallest(xmean,samples),largest(xmean,samples),"<fluoro>","#",sim);
      
      if(AutoConditioningLevelQ) {
        sim.io <<"guessing optimal conditioning level..."<<Endl;
        GlobalConditioningLevel = Magic_GuessConditioningLevel(xdatadouble,size,samples,sim);
        sim.io <<" -> conditioning level is: "<<GlobalConditioningLevel<<Endl;
        sim.io <<" -> done."<<Endl;
      } else {
        if(GlobalConditioningLevel >= 0.0 && RelativeGlobalConditioningLevelQ) {
          sim.io <<"using relative conditioning level, determining absolute value..."<<Endl;
          double xmean_min = smallest(xmean, samples);
          double xmean_max = largest(xmean, samples);
          // rescale conditioning level
          GlobalConditioningLevel = (xmean_max - xmean_min) * GlobalConditioningLevel + xmean_min;
          sim.io <<" -> conditioning level is now: "<<GlobalConditioningLevel<<Endl;
        }
      }

      free_time_series_memory(xmean);
      
      if((globalbins>1)&&(!GenerateGlobalFromFilteredDataQ)) {
        sim.io <<"generating global signal and collecting local signals..."<<Endl;
        xglobaldouble = generate_mean_time_series(xdatadouble, size, samples);
        xdatadoubleglue = generate_conditioned_time_series_by_glueing(xdatadouble, size, xglobaldouble, StartSampleIndex, EndSampleIndex, GlobalConditioningLevel, &AvailableSamples, sim);
        sim.io <<" -> done."<<Endl;
      }
      
      if(HighPassFilterQ) {
        sim.io <<"applying high-pass filter to time series..."<<Endl;
        apply_high_pass_filter_to_time_series(xdatadouble, size, samples);
        if(globalbins>1)
          apply_high_pass_filter_to_time_series(xdatadoubleglue, size, AvailableSamples);
        sim.io <<" -> done."<<Endl;
      }

      if((globalbins>1)&&GenerateGlobalFromFilteredDataQ) {
        sim.io <<"generating global signal and collecting local signals..."<<Endl;
        xglobaldouble = generate_mean_time_series(xdatadouble, size, samples);
        sim.io <<"... step 2"<<Endl;
        xdatadoubleglue = generate_conditioned_time_series_by_glueing(xdatadouble, size, xglobaldouble, StartSampleIndex, EndSampleIndex, GlobalConditioningLevel, &AvailableSamples, sim);
        sim.io <<" -> done."<<Endl;
      }

      sim.io <<" -> done."<<Endl;
    }
    catch(...) {
      sim.io <<"Error: could not reserve enough memory!"<<Endl;
      if(!ContinueOnErrorQ) exit(1);
      else
      {
        sim.io <<"Error handling: ContinueOnErrorQ flag set, proceeding..."<<Endl;
        skip_the_rest = true;
      }
    }
    
    // if (AvailableSamples < 1) {
    //       sim.io <<"Warning: No samples available, skipping evaluation!"<<Endl;
    //       skip_the_rest = true;
    // }
    
    if (!skip_the_rest) {
    // main loop:
    sim.io <<"set-up: "<<size<<" nodes, ";
    sim.io <<EndSampleIndex-StartSampleIndex+1<<" out of "<<samples<<" samples, ";
  
    time(&start);
    sim.io <<"start: "<<ctime(&start)<<Endl;
#ifdef SHOW_DETAILED_PROGRESS
    sim.io <<"running ";
#else
    sim.io <<"running..."<<Endl;
    bool status_already_displayed = false;
#endif

    for(int ii=0; ii<size; ii++)
    {
#ifdef SHOW_DETAILED_PROGRESS
      status(ii,REPORTS,size);
#else
      time(&middle);
      if ((!status_already_displayed)&&((ii>=size/3)||((middle-start>30.)&&(ii>0))))
      { 
        sim.io <<" (after "<<ii<<" nodes: elapsed "<<sec2string(difftime(middle,start)) \
          <<", ETA "<<ETAstring(ii,size,difftime(middle,start))<<")"<<Endl;
        status_already_displayed = true;
      }
#endif
      for(int jj=0; jj<size; jj++)
      {
        if (ii != jj)
        {
#ifndef SEPARATED_OUTPUT
          if(globalbins>1)
            xresult[jj][ii] = CrossCorrelation(xdatadoubleglue[ii], xdatadoubleglue[jj]);
          else xresult[jj][ii] = CrossCorrelation(xdatadouble[ii], xdatadouble[jj]);
#else
          if(globalbins>1)
            CrossCorrelationSeparated(xdatadoubleglue[ii], xdatadoubleglue[jj], ii, jj);
          CrossCorrelationSeparated(xdatadouble[ii], xdatadouble[jj], ii, jj);
#endif
        }
        // else xresult[ii][jj] = 0.0;
      }
    }
#ifndef SHOW_DETAILED_PROGRESS
    sim.io <<" -> done."<<Endl;
#endif

    time(&end);
    sim.io <<"end: "<<ctime(&end)<<Endl;
    sim.io <<"runtime: "<<sec2string(difftime(end,start))<<Endl;

    // cout <<"TE terms: "<<terms_sum<<", of those zero: "<<terms_zero<<" ("<<int(double(terms_zero)*100/terms_sum)<<"%)"<<endl;
    
  }
  };
  
  void finalize(Sim& sim)
  {
    if(!skip_the_rest) {
#ifdef SEPARATED_OUTPUT
      write_multidim_result(xresult,globalbins,size,outputfile_results_name,sim);
#else
      write_result(xresult,size,outputfile_results_name,sim);
#endif
      write_parameters();

      gsl_rng_free(GSLrandom);
    }

    try {
#ifdef SEPARATED_OUTPUT
    for (int x=0; x<size; x++)
    {
      for (int x2=0; x2<size; x2++)
        delete[] xresult[x][x2];
      delete[] xresult[x];
    }
#endif
    delete[] xresult;
    
    // if (AvailableSamples != NULL) delete[] AvailableSamples;

    // free_time_series_memory(xdatadouble,size); earlier...
    if (xdatadouble!=NULL) free_time_series_memory(xdatadouble,size);
    if (xdatadoubleglue!=NULL) free_time_series_memory(xdatadoubleglue,size);
    if (xglobaldouble!=NULL) free_time_series_memory(xglobaldouble);
    }
    catch(...) {};
    
    sim.io <<"End of Kernel (iteration="<<(sim.iteration())<<")"<<Endl;
  };
  

#ifdef SEPARATED_OUTPUT
  void CrossCorrelationSeparated(double *arrayI, double *arrayJ, int I, int J)
#else
  double CrossCorrelation(double *arrayI, double *arrayJ)
#endif
  {
    // We are looking at the information flow of array1 ("J") -> array2 ("I")
    double result = 0.0;
    double temp;
    unsigned long samples_to_use_here = (unsigned long)(samples);
    unsigned long const JShift = (unsigned long const)InstantFeedbackTermQ;
    if (GlobalConditioningLevel>=0.) samples_to_use_here = AvailableSamples;
    assert(AvailableSamples<=EndSampleIndex-StartSampleIndex+1);

    if(samples_to_use_here>CROSSCORRELATION_MAX_LAG+2) {
      // include lag of zero depending on InstantFeedbackTermQ flag
      for(int lag=1-JShift;lag<=CROSSCORRELATION_MAX_LAG; lag++) {
        temp = gsl_stats_correlation(&arrayI[lag],1,&arrayJ[0],1,samples_to_use_here-lag);
        if (temp>result) result = temp;
      }
    }

    return result;
  };

  void write_parameters()
  {
    char* name = new char[outputfile_pars_name.length()+1];
    strcpy(name,outputfile_pars_name.c_str());
    ofstream fileout1(name);
    delete[] name;
    if (fileout1 == NULL)
    {
      cerr <<endl<<"error: cannot open output file!"<<endl;
      exit(1);
    }

    fileout1.precision(6);
    fileout1 <<"{";
    fileout1 <<"executable->xcorrsim";
    fileout1 <<", iteration->"<<iteration;
    time(&end);
    fileout1 <<", ExecutionTime->\""<<sec2string(difftime(end,start))<<"\"";

    fileout1 <<", size->"<<size;
    // fileout1 <<", bins->"<<bins;
    fileout1 <<", globalbins->"<<globalbins;
    fileout1 <<", appliedscaling->"<<input_scaling;
    fileout1 <<", samples->"<<samples;
    fileout1 <<", StartSampleIndex->"<<StartSampleIndex;
    fileout1 <<", EndSampleIndex->"<<EndSampleIndex;
    fileout1 <<", EqualSampleNumberQ->"<<bool2textMX(EqualSampleNumberQ);
    fileout1 <<", MaxSampleNumberPerBin->"<<MaxSampleNumberPerBin;
    fileout1 <<", AvailableSamples->{"<<AvailableSamples<<"}";

    fileout1 <<", cutoff->"<<cutoff;
    fileout1 <<", noise->"<<std_noise;
    fileout1 <<", tauF->"<<tauF;
    fileout1 <<", tauCa->"<<tauCa;
    fileout1 <<", DeltaCalciumOnAP->"<<DeltaCalciumOnAP;
    fileout1 <<", OverrideRescalingQ->"<<bool2textMX(OverrideRescalingQ);
    fileout1 <<", HighPassFilterQ->"<<bool2textMX(HighPassFilterQ);
    fileout1 <<", InstantFeedbackTermQ->"<<bool2textMX(InstantFeedbackTermQ);
// #ifdef ENABLE_ADAPTIVE_BINNING_AT_COMPILE_TIME
//    fileout1 <<", AdaptiveBinningQ->"<<bool2textMX(AdaptiveBinning)Q;
// #endif
    fileout1 <<", ContinueOnErrorQ->"<<bool2textMX(ContinueOnErrorQ);
    fileout1 <<", saturation->"<<fluorescence_saturation;
    fileout1 <<", IncludeGlobalSignalQ->"<<bool2textMX(IncludeGlobalSignalQ);
    fileout1 <<", GenerateGlobalFromFilteredDataQ->"<<bool2textMX(GenerateGlobalFromFilteredDataQ);
    fileout1 <<", GlobalConditioningLevel->"<<GlobalConditioningLevel;
    // fileout1 <<", TargetMarkovOrder->"<<TargetMarkovOrder;
    // fileout1 <<", SourceMarkovOrder->"<<SourceMarkovOrder;

    // fileout1 <<", AutoBinNumberQ->"<<bool2textMX(AutoBinNumberQ);
    fileout1 <<", AutoConditioningLevelQ->"<<bool2textMX(AutoConditioningLevelQ);


    fileout1 <<", inputfile->\""<<inputfile_name<<"\"";
    fileout1 <<", outputfile->\""<<outputfile_results_name<<"\"";
    fileout1 <<", spikeindexfile->\""<<spikeindexfile_name<<"\"";
    fileout1 <<", spiketimesfile->\""<<spiketimesfile_name<<"\"";
    fileout1 <<", FluorescenceModel->\""<<FluorescenceModel<<"\"";
    // parameters for light scattering
    fileout1 <<", YAMLfile->\""<<YAMLfilename<<"\"";
    fileout1 <<", SigmaScatter->"<<SigmaScatter;
    fileout1 <<", AmplitudeScatter->"<<AmplitudeScatter;
    fileout1 <<"}"<<endl;

    fileout1.close();
  };
  
};
