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

// Calculate the Mutual Information between a number of time series.

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <ctime>
#include <cstring>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_vector.h>

#include "../olav.h"
#include <sim_main.h>
#include "../te-datainit.h"
#include "../multidimarray.h"


#define REPORTS 25

#define OUTPUTNUMBER_PRECISION 15
#define SEPARATED_OUTPUT

#define MUTUALINFORMATION_MAX_LAG 4

// #define GSL_RANDOM_NUMBER_GENERATOR gsl_rng_default
#define GSL_RANDOM_NUMBER_GENERATOR gsl_rng_ranlxs2

#define COUNTARRAY_IPAST_GPAST 1
#define COUNTARRAY_JPAST_GPAST 2
#define COUNTARRAY_IPAST_JPAST_GPAST 3

using namespace std;

typedef unsigned char rawdata;

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
  unsigned int bins, globalbins;
  // unsigned long mag der SimKernel irgendwie nicht?
  long samples;
  long StartSampleIndex, EndSampleIndex;
  bool EqualSampleNumberQ;
  long MaxSampleNumberPerBin;
  unsigned long * AvailableSamples;
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

  bool IncludeGlobalSignalQ;
  bool GenerateGlobalFromFilteredDataQ;
  double GlobalConditioningLevel;
  int SourceMarkovOrder, TargetMarkovOrder;

  bool ContinueOnErrorQ;
  bool skip_the_rest;

  bool AutoBinNumberQ;
  bool AutoConditioningLevelQ;

  // Orlandi: Adding option for predefined binning limits
  std::vector<double> binEdges;
  bool binEdgesSet;

  gsl_rng* GSLrandom;

  unsigned long Tspace, Sspace;
  MultiDimArrayLong* F_Ipast_Gpast;
  MultiDimArrayLong* F_Jpast_Gpast;
  MultiDimArrayLong* F_Ipast_Jpast_Gpast;

  gsl_vector_int * vec_Full;
  gsl_vector_int * vec_Full_Bins;
  
  gsl_vector_int_view vec_Ipast;
  gsl_vector_int_view vec_Jpast;
  // here the conditioning signal is fixed to order 1
  gsl_vector_int_view vec_Gpast;
  gsl_vector_int* gsl_access;
  
  rawdata **xdata;
  rawdata *xglobal;
#ifndef SEPARATED_OUTPUT
  double **xresult;
  long double Hxx, Hxxy;
#else
  double ***xresult;
  long double *MItemp;
  long double *MIresult;
#endif

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
    bins = 0;
    sim.get("AutoBinNumberQ",AutoBinNumberQ,false);
    if(!AutoBinNumberQ) sim.get("bins",bins);
    sim.get("globalbins",globalbins,1);
    if(globalbins<1) globalbins=1;
    sim.get("samples",samples);
    sim.get("StartSampleIndex",StartSampleIndex,1);
    sim.get("EndSampleIndex",EndSampleIndex,samples-1);
    sim.get("EqualSampleNumberQ",EqualSampleNumberQ,false);
    sim.get("MaxSampleNumberPerBin",MaxSampleNumberPerBin,-1);

    // Orlandi: Adding option for predefined binning limits
    sim.get("binEdges", binEdges, 0);
    if(binEdges.size() > 1)
    {
      // Check that the bin edges are increasing in size
      for(std::vector<double>::iterator iteratorBinEdges = binEdges.begin()+1; iteratorBinEdges != binEdges.end(); iteratorBinEdges++)
        assert(*iteratorBinEdges-*(iteratorBinEdges-1) > 0.);
      assert(binEdges.size()-1 == bins);
      binEdgesSet = true;
    }
    else
      binEdgesSet = false;

    sim.get("noise",std_noise,-1.);
    sim.get("appliedscaling",input_scaling,1.);
    sim.get("cutoff",cutoff,-1.);
    sim.get("tauF",tauF,-1.);
    sim.get("OverrideRescalingQ",OverrideRescalingQ,false);
    sim.get("HighPassFilterQ",HighPassFilterQ,false);
    sim.get("InstantFeedbackTermQ",InstantFeedbackTermQ,false);

    sim.get("saturation",fluorescence_saturation,-1.);
    sim.get("IncludeGlobalSignalQ",IncludeGlobalSignalQ,false);
    assert(IncludeGlobalSignalQ ^ globalbins==1);
    sim.get("GenerateGlobalFromFilteredDataQ",GenerateGlobalFromFilteredDataQ,false);
    sim.get("AutoConditioningLevelQ",AutoConditioningLevelQ,false);
    if(!AutoConditioningLevelQ)
    {
      sim.get("GlobalConditioningLevel",GlobalConditioningLevel,-1.);
      if (GlobalConditioningLevel>0) assert(globalbins==2); // for now, to keep it simple
    }
    else GlobalConditioningLevel = -1.;
    
    sim.get("SourceMarkovOrder",SourceMarkovOrder,1);
    assert(SourceMarkovOrder>0);
    sim.get("TargetMarkovOrder",TargetMarkovOrder,1);
    assert(TargetMarkovOrder>0);
    
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
    
    AvailableSamples = NULL;
    xdata = NULL;
    xglobal = NULL;
    xresult = NULL;
    
#ifdef SEPARATED_OUTPUT
    MItemp = NULL;
    MIresult = NULL;
#endif
  };

  void execute(Sim& sim)
  {
    sim.io <<"------ transferentropy-sim:extended-multidim ------ olav, Wed 7 Sep 2011 ------"<<Endl;

    sim.io <<"output file: "<<outputfile_results_name<<Endl;
    // Gespeichert wird später - hier nur Test, ob das Zielverzeichnis existiert
    write_parameters();
    skip_the_rest = false;

    sim.io <<"allocating memory..."<<Endl;
    try {
#ifndef SEPARATED_OUTPUT
      xresult = new double*[size];
#else
      xresult = new double**[size];
#endif
      for(int i=0; i<size; i++)
      {
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
    
#ifdef SEPARATED_OUTPUT
      MItemp = new long double[globalbins];
      MIresult = new long double[globalbins];
#endif
      AvailableSamples = new unsigned long[globalbins];
      
      // hack of medium ugliness to make it work without global signal
      if(globalbins<=1)
      {
        // sim.io <<"debug: xglobal hack."<<Endl;
        xglobal = new rawdata[samples];
        memset(xglobal, 0, samples*sizeof(rawdata));
        AvailableSamples[0] = EndSampleIndex-StartSampleIndex+1;
      }
      sim.io <<" -> done."<<Endl;     
      
      double** xdatadouble = NULL;
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

      // sim.io <<"histogram of averaged signal:"<<Endl;
      // double* xmean = generate_mean_time_series(xdatadouble,size,samples);
      // PlotLogHistogramInASCII(xmean,samples,smallest(xmean,samples),largest(xmean,samples),"<fluoro>","#",sim);
      // free_time_series_memory(xmean);
      // cout <<"DEBUG: subset of first node: ";
      // display_subset(xdatadouble[0]);
      if(AutoConditioningLevelQ) {
        sim.io <<"guessing optimal conditioning level..."<<Endl;
        GlobalConditioningLevel = Magic_GuessConditioningLevel(xdatadouble,size,samples,sim);
        sim.io <<" -> conditioning level is: "<<GlobalConditioningLevel<<Endl;
        sim.io <<" -> done."<<Endl;
      }      
      
      if((globalbins>1)&&(!GenerateGlobalFromFilteredDataQ)) {
        sim.io <<"generating discretized global signal..."<<Endl;
        xglobal = generate_discretized_global_time_series(xdatadouble, size, samples, globalbins, GlobalConditioningLevel, AvailableSamples, StartSampleIndex, EndSampleIndex, EqualSampleNumberQ, MaxSampleNumberPerBin, sim);
        sim.io <<" -> done."<<Endl;
      }
      
      if(HighPassFilterQ) {
        sim.io <<"applying high-pass filter to time series..."<<Endl;
        apply_high_pass_filter_to_time_series(xdatadouble, size, samples);
        sim.io <<" -> done."<<Endl;
      }

      if(AutoBinNumberQ) {
        sim.io <<"guessing optimal bin number..."<<Endl;
        bins = Magic_GuessBinNumber(xdatadouble,size,samples);
        sim.io <<" -> number of bins is: "<<bins<<Endl;
        sim.io <<" -> done."<<Endl;
      }
      
      // Now we know the number of local bins to use, so we can reserve the discretized memory:
      // This is overall iterator that will be mapped onto array indices later:
      vec_Full = NULL;
      vec_Full_Bins = NULL;
      vec_Full = gsl_vector_int_alloc(TargetMarkovOrder+SourceMarkovOrder+1);
      gsl_vector_int_set_zero(vec_Full);
      vec_Full_Bins = gsl_vector_int_alloc(TargetMarkovOrder+SourceMarkovOrder+1);
      gsl_vector_int_set_all(vec_Full_Bins,bins);
      gsl_vector_int_set(vec_Full_Bins,TargetMarkovOrder+SourceMarkovOrder,globalbins);
      gsl_access = gsl_vector_int_alloc(TargetMarkovOrder+SourceMarkovOrder+1);
    
      // Initialize views to have better access to the full iterator elements:
      vec_Ipast = gsl_vector_int_subvector(vec_Full,0,TargetMarkovOrder);
      vec_Jpast = gsl_vector_int_subvector(vec_Full,TargetMarkovOrder,SourceMarkovOrder);
      vec_Gpast = gsl_vector_int_subvector(vec_Full,TargetMarkovOrder+SourceMarkovOrder,1);
      
      // here we assume equal binning for source and target terms
      // ------------------ IndexMultipliers_Ipast_Gpast:
      gsl_vector_int* BinsPerDim = gsl_vector_int_alloc(TargetMarkovOrder+1);
      for (int i=0; i<TargetMarkovOrder; i++)
        gsl_vector_int_set(BinsPerDim,i,bins);
      gsl_vector_int_set(BinsPerDim,TargetMarkovOrder,globalbins);
      F_Ipast_Gpast = new MultiDimArrayLong(BinsPerDim);
      gsl_vector_int_free(BinsPerDim);

      // ------------------ IndexMultipliers_Jpast_Gpast:
      BinsPerDim = gsl_vector_int_alloc(SourceMarkovOrder+1);
      for (int i=0; i<SourceMarkovOrder; i++)
        gsl_vector_int_set(BinsPerDim,i,bins);
      gsl_vector_int_set(BinsPerDim,SourceMarkovOrder,globalbins);
      F_Jpast_Gpast = new MultiDimArrayLong(BinsPerDim);
      gsl_vector_int_free(BinsPerDim);
    
      // ------------------ IndexMultipliers_Ipast_Jpast_Gpast:
      BinsPerDim = gsl_vector_int_alloc(TargetMarkovOrder+SourceMarkovOrder+1);
      for (int i=0; i<TargetMarkovOrder+SourceMarkovOrder; i++)
        gsl_vector_int_set(BinsPerDim,i,bins);
      gsl_vector_int_set(BinsPerDim,TargetMarkovOrder+SourceMarkovOrder,globalbins);
      F_Ipast_Jpast_Gpast = new MultiDimArrayLong(BinsPerDim);
      gsl_vector_int_free(BinsPerDim);
    
           
      if((globalbins>1)&&(GenerateGlobalFromFilteredDataQ)) {
        sim.io <<"generating discretized global signal..."<<Endl;
        xglobal = generate_discretized_global_time_series(xdatadouble, size, samples, globalbins, GlobalConditioningLevel, AvailableSamples, StartSampleIndex, EndSampleIndex, EqualSampleNumberQ, MaxSampleNumberPerBin, sim);
        sim.io <<" -> done."<<Endl;
      }

      sim.io <<"discretizing local time series..."<<Endl;
      if(!binEdgesSet)
        xdata = generate_discretized_version_of_time_series(xdatadouble, size, samples, bins);
      else
        xdata = generate_discretized_version_of_time_series(xdatadouble, size, samples, binEdges);
      // and, since double version of time series is not used any more...
      if(xdatadouble!=NULL) free_time_series_memory(xdatadouble, size);
      sim.io <<" -> done."<<Endl;

      // cout <<"DEBUG: subset of discretized global signal: ";
      // display_subset(xglobal);
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
    
    if (!skip_the_rest) {
  
    // main loop:
    sim.io <<"set-up: "<<size<<" nodes, ";
    sim.io <<EndSampleIndex-StartSampleIndex+1<<" out of "<<samples<<" samples, ";
    sim.io <<bins<<" bins, "<<globalbins<<" globalbins"<<Endl;
    sim.io <<"set-up: Markov order of source/target/conditioning: "<<SourceMarkovOrder<<"/"<<TargetMarkovOrder<<"/1"<<Endl;
#ifdef SEPARATED_OUTPUT
    sim.io <<"set-up: separated output (globalbin)"<<Endl;
#endif
  
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
          xresult[jj][ii] = MutualInformation(xdata[ii], xdata[jj]);
#else
          MutualInformationSeparated(xdata[ii], xdata[jj], ii, jj);
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

  }
  };
  
  void finalize(Sim& sim)
  {
    if(!skip_the_rest) {
#ifdef SEPARATED_OUTPUT
      write_multidim_result(xresult,globalbins);
#else
      write_result(xresult);
#endif
      write_parameters();

      // free allocated memory
      gsl_rng_free(GSLrandom);
      gsl_vector_int_free(vec_Full);
      gsl_vector_int_free(vec_Full_Bins);
      gsl_vector_int_free(gsl_access);
    }

    try {
#ifdef SEPARATED_OUTPUT
    for (int x=0; x<size; x++)
    {
      for (int x2=0; x2<size; x2++)
        delete[] xresult[x][x2];
      delete[] xresult[x];
    }
    if(MItemp != NULL) delete[] MItemp;
    if(MIresult != NULL) delete[] MIresult;
#endif
    delete[] xresult;
    
    if (AvailableSamples != NULL) delete[] AvailableSamples;

    delete F_Ipast_Gpast;
    delete F_Jpast_Gpast;
    delete F_Ipast_Jpast_Gpast;

    if(xdata != NULL) free_time_series_memory(xdata,size);
    if(xglobal != NULL) free_time_series_memory(xglobal);
    }
    catch(...) {};
    
    sim.io <<"End of Kernel (iteration="<<(sim.iteration())<<")"<<Endl;
  };
  

#ifdef SEPARATED_OUTPUT
  void MutualInformationSeparated(rawdata *arrayI, rawdata *arrayJ, int I, int J)
#else
  double MutualInformation(rawdata *arrayI, rawdata *arrayJ)
#endif
  {  
    bool runningIt;
    unsigned long ig, jg, ijg;
    long double igd, jgd, ijgd;
#ifndef SEPARATED_OUTPUT
    double MIresult = 0.0;
    double MItemp;
#else
    memset(MIresult,0,globalbins*sizeof(long double));
#endif

    for(long lag=0;lag<=MUTUALINFORMATION_MAX_LAG; lag++) // start lag loop
    {
#ifdef SEPARATED_OUTPUT
      memset(MItemp,0,globalbins*sizeof(long double));
#else
      tempresult = 0.;
#endif
      F_Ipast_Gpast->clear();
      F_Jpast_Gpast->clear();
      F_Ipast_Jpast_Gpast->clear();

      // unsigned long const JShift = (unsigned long const)InstantFeedbackTermQ;
      assert(StartSampleIndex >= max(TargetMarkovOrder,SourceMarkovOrder));
  
      // extract probabilities (actually number of occurrence)
      for (unsigned long t=max(StartSampleIndex,lag); t<=EndSampleIndex; t++)
      {
        // prepare the index vector vec_Full via the vector views
        for (int i=0; i<TargetMarkovOrder; i++)
          gsl_vector_int_set(&vec_Ipast.vector,i,arrayI[t-i]);
        
        for (int i=0; i<SourceMarkovOrder; i++)
          gsl_vector_int_set(&vec_Jpast.vector,i,arrayJ[t-i-lag]);
        
        gsl_vector_int_set(&vec_Gpast.vector,0,xglobal[t]);
      
        // add counts to arrays
        if (xglobal[t]<globalbins) { // for EqualSampleNumberQ flag
          set_up_access_vector(COUNTARRAY_IPAST_GPAST);
          F_Ipast_Gpast->inc(gsl_access);
          set_up_access_vector(COUNTARRAY_JPAST_GPAST);
          F_Jpast_Gpast->inc(gsl_access);
          set_up_access_vector(COUNTARRAY_IPAST_JPAST_GPAST);
          F_Ipast_Jpast_Gpast->inc(gsl_access);
        }
      }
    
      gsl_vector_int_set_zero(vec_Full);
      runningIt = true;
      while (runningIt)
      {
        // SimplePrintGSLVector(vec_Full);
        set_up_access_vector(COUNTARRAY_IPAST_JPAST_GPAST);
        ijg = F_Ipast_Jpast_Gpast->get(gsl_access);
        ijgd = (long double)ijg;
        if (ijg!=0)
        {
          set_up_access_vector(COUNTARRAY_IPAST_GPAST);
          ig = F_Ipast_Gpast->get(gsl_access);
          igd = (long double)ig;
          set_up_access_vector(COUNTARRAY_JPAST_GPAST);
          jg = F_Jpast_Gpast->get(gsl_access);
          jgd = (long double)jg;
        
#ifndef SEPARATED_OUTPUT
          MItemp += 
#else
          MItemp[gsl_vector_int_get(&vec_Gpast.vector,0)] =
#endif
            ijgd/AvailableSamples[gsl_vector_int_get(&vec_Gpast.vector,0)] * \
            log(ijgd*AvailableSamples[gsl_vector_int_get(&vec_Gpast.vector,0)]/(igd*jgd));
      
        }
        runningIt = OneStepAhead_FullIterator();
      }
      
      // see if the MI values are greater for the current lag, and save them for later if so
#ifndef SEPARATED_OUTPUT
      if(MItemp > MIresult) MIresult = MItemp;
#else
      for (rawdata g=0; g<globalbins; g++) {
        if(MItemp[g] > MIresult[g]) MIresult[g] = MItemp[g];
      }
#endif
    } // end lag loop
    
#ifdef SEPARATED_OUTPUT
    for (rawdata g=0; g<globalbins; g++) {
      xresult[J][I][g] = double(MIresult[g]/log(2.));
    }
#else
    return double((Hxx - Hxxy)/log(2));
#endif
  };

  std::string bool2textMX(bool value)
  {
    if (value) return "True";
    else return "False";
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
    fileout1 <<"executable->mi";
    fileout1 <<", iteration->"<<iteration;
    time(&end);
    fileout1 <<", ExecutionTime->"<<sec2string(difftime(end,start));
    
    fileout1 <<", size->"<<size;
    fileout1 <<", bins->"<<bins;
    fileout1 <<", globalbins->"<<globalbins;
    fileout1 <<", appliedscaling->"<<input_scaling;
    fileout1 <<", samples->"<<samples;
    fileout1 <<", StartSampleIndex->"<<StartSampleIndex;
    fileout1 <<", EndSampleIndex->"<<EndSampleIndex;
    fileout1 <<", EqualSampleNumberQ->"<<bool2textMX(EqualSampleNumberQ);
    fileout1 <<", MaxSampleNumberPerBin->"<<MaxSampleNumberPerBin;
    fileout1 <<", AvailableSamples->{";
    for (int i=0; i<globalbins; i++)
    {
      if (i>0) fileout1 <<",";
      if (AvailableSamples == NULL) fileout1 <<"?";
      else fileout1 <<AvailableSamples[i];
    }
    fileout1 <<"}";

    fileout1 <<", cutoff->"<<cutoff;
    fileout1 <<", noise->"<<std_noise;
    fileout1 <<", tauF->"<<tauF;
    fileout1 <<", tauCa->"<<tauCa;
    fileout1 <<", DeltaCalciumOnAP->"<<DeltaCalciumOnAP;
    fileout1 <<", OverrideRescalingQ->"<<bool2textMX(OverrideRescalingQ);
    fileout1 <<", HighPassFilterQ->"<<bool2textMX(HighPassFilterQ);
    fileout1 <<", InstantFeedbackTermQ->"<<bool2textMX(InstantFeedbackTermQ);

    fileout1 <<", ContinueOnErrorQ->"<<bool2textMX(ContinueOnErrorQ);
    fileout1 <<", saturation->"<<fluorescence_saturation;
    fileout1 <<", IncludeGlobalSignalQ->"<<bool2textMX(IncludeGlobalSignalQ);
    fileout1 <<", GenerateGlobalFromFilteredDataQ->"<<bool2textMX(GenerateGlobalFromFilteredDataQ);
    fileout1 <<", GlobalConditioningLevel->"<<GlobalConditioningLevel;
    fileout1 <<", TargetMarkovOrder->"<<TargetMarkovOrder;
    fileout1 <<", SourceMarkovOrder->"<<SourceMarkovOrder;
    
    fileout1 <<", AutoBinNumberQ->"<<bool2textMX(AutoBinNumberQ);
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

  void write_result(double **array)
  {
    char* name = new char[outputfile_results_name.length()+1];
    strcpy(name,outputfile_results_name.c_str());
    ofstream fileout1(name);
    delete[] name;
    if (fileout1 == NULL)
    {
      cerr <<endl<<"error: cannot open output file!"<<endl;
      exit(1);
    }   

    fileout1.precision(OUTPUTNUMBER_PRECISION);
    fileout1 <<fixed;
    fileout1 <<"{";
    for(int j=0; j<size; j++)
    {
      if(j>0) fileout1<<",";
      fileout1 <<"{";
      for(unsigned long i=0; i<size; i++)
      {
        if (i>0) fileout1<<",";
        fileout1 <<(double)array[j][i];
      }
      fileout1 <<"}"<<endl;
    }
    fileout1 <<"}"<<endl;

    cout <<"Transfer entropy matrix saved."<<endl;
  };

  void write_multidim_result(double ***array, unsigned int dimens)
  {
    char* name = new char[outputfile_results_name.length()+1];
    strcpy(name,outputfile_results_name.c_str());
    ofstream fileout1(name);
    delete[] name;
    if (fileout1 == NULL)
    {
      cerr <<endl<<"error: cannot open output file!"<<endl;
      exit(1);
    }   

    fileout1.precision(OUTPUTNUMBER_PRECISION);
    fileout1 <<fixed;
    fileout1 <<"{";
    for(unsigned long j=0; j<size; j++)
    {
      if(j>0) fileout1<<",";
      fileout1 <<"{";
      for(unsigned long i=0; i<size; i++)
      {
        if (i>0) fileout1<<",";
        fileout1 <<"{";
        for(int k=0; k<dimens; k++)
        {
          if (k>0) fileout1<<",";
          fileout1 <<array[j][i][k];
        }
        fileout1 <<"}";
      }
      fileout1 <<"}"<<endl;
    }
    fileout1 <<"}"<<endl;

    cout <<"Transfer entropy matrix saved."<<endl;
  };

  bool OneStepAhead_FullIterator()
  {
    for(int i=0; i<=TargetMarkovOrder+SourceMarkovOrder+1; i++)
    {
      if (i==TargetMarkovOrder+SourceMarkovOrder+1) return false; // if we have reached the "maximum" value
      gsl_vector_int_set(vec_Full,i,gsl_vector_int_get(vec_Full,i)+1);
      
      if (gsl_vector_int_get(vec_Full,i) >= gsl_vector_int_get(vec_Full_Bins,i))
      {
        gsl_vector_int_set(vec_Full,i,0);
        if((i==TargetMarkovOrder+SourceMarkovOrder)&&(GlobalConditioningLevel>0.0))
          return false; // if we have reached the effective maximum, because we don't want to go through more
      }
      else break;
    }
    return true;
  };
  
  void SimplePrintGSLVector(gsl_vector_int* vec, Sim& sim)
  {
    SimplePrintGSLVector(vec,true,sim);
  };
  void SimplePrintGSLVector(gsl_vector_int* vec, bool newline, Sim& sim)
  {
    for(int i=0; i<vec->size; i++)
      sim.io <<gsl_vector_int_get(vec,i)<<" ";
    if (newline) sim.io <<Endl;
  };

  void SimplePrintFullIterator()
  {
    SimplePrintFullIterator(true);
  };
  void SimplePrintFullIterator(bool newline)
  {
    // cout <<"Inow: "<<gsl_vector_int_get(&vec_Inow.vector,0);
    cout <<" Ipast: ";
    for (int i=0; i<TargetMarkovOrder; i++)
      cout <<gsl_vector_int_get(&vec_Ipast.vector,i)<<" ";
    cout <<"Jpast: ";
    for (int i=0; i<SourceMarkovOrder; i++)
      cout <<gsl_vector_int_get(&vec_Jpast.vector,i)<<" ";
    cout <<"Gpast: "<<gsl_vector_int_get(&vec_Gpast.vector,0);
    if (newline) cout <<endl;
  };
  
  // The following assumes that vec_Full has been set already (!)
  void set_up_access_vector(int arraycode)
  {
    switch (arraycode)
    {
      case COUNTARRAY_IPAST_GPAST:
        for (int i=0; i<TargetMarkovOrder; i++)
          gsl_vector_int_set(gsl_access,i,gsl_vector_int_get(&vec_Ipast.vector,i));
        gsl_vector_int_set(gsl_access,TargetMarkovOrder,gsl_vector_int_get(&vec_Gpast.vector,0));
        break;

      case COUNTARRAY_JPAST_GPAST:
        for (int i=0; i<SourceMarkovOrder; i++)
          gsl_vector_int_set(gsl_access,i,gsl_vector_int_get(&vec_Jpast.vector,i));
        gsl_vector_int_set(gsl_access,SourceMarkovOrder,gsl_vector_int_get(&vec_Gpast.vector,0));
        break;

      case COUNTARRAY_IPAST_JPAST_GPAST:
        for (int i=0; i<TargetMarkovOrder; i++)
          gsl_vector_int_set(gsl_access,i,gsl_vector_int_get(&vec_Ipast.vector,i));
        for (int i=0; i<SourceMarkovOrder; i++)
          gsl_vector_int_set(gsl_access,TargetMarkovOrder+i,gsl_vector_int_get(&vec_Jpast.vector,i));
        gsl_vector_int_set(gsl_access,SourceMarkovOrder+TargetMarkovOrder,gsl_vector_int_get(&vec_Gpast.vector,0));
        break;

      default:
        cout <<endl<<"GetCounterArray: error, invalid array code."<<endl; exit(1);
    }
  }
};
