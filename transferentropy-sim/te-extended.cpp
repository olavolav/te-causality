// calculate the transfer entropy between a numer of time series
// this is the extension to arbitrary Markov order of the source and target term
// created by olav, Di 12 Okt 2010 11:34:30 CEST

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

#include "../../Simulationen/olav.h"
#include "../../../Sonstiges/SimKernel/sim_main.h"
#include "../te-datainit.h"

// #ifndef INLINE
// #define INLINE extern inline
// #endif

#define REPORTS 25
#undef SHOW_DETAILED_PROGRESS

#define OUTPUTNUMBER_PRECISION 15
#define SEPARATED_OUTPUT

#undef ENABLE_ADAPTIVE_BINNING_AT_COMPILE_TIME

#define COUNTARRAY_IPAST_GPAST 1
#define COUNTARRAY_INOW_IPAST_GPAST 2
#define COUNTARRAY_IPAST_JPAST_GPAST 3
#define COUNTARRAY_INOW_IPAST_JPAST_GPAST 4

using namespace std;

typedef unsigned char rawdata;

time_t start, now;

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
	unsigned int rawdatabins;
	// unsigned long mag der SimKernel irgendwie nicht?
	long samples;
	long StartSampleIndex, EndSampleIndex;
	bool EqualSampleNumberQ;
	long MaxSampleNumberPerBin;
	unsigned long * AvailableSamples;
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
  
	bool OverrideRescalingQ; // if, for example the input are pure spike data (integers)
	bool HighPassFilterQ; // actually, this transforms the signal into the difference signal
	bool InstantFeedbackTermQ;
#ifdef ENABLE_ADAPTIVE_BINNING_AT_COMPILE_TIME
	bool AdaptiveBinningQ; // rubbish
#endif
	bool IncludeGlobalSignalQ;
	bool GenerateGlobalFromFilteredDataQ;
	double GlobalConditioningLevel;
	int SourceMarkovOrder, TargetMarkovOrder;
	// speed test:
	bool SingleSourceMarkovOrder, SingleTargetMarkovOrder;
	
	bool ContinueOnErrorQ;
	bool skip_the_rest;
	
  bool AutoBinNumberQ;
  bool AutoConditioningLevelQ;

	// using the new smart gsl_vector things, the F_arrays are all 1D internally
	unsigned long * F_Ipast_Gpast;
	unsigned long * F_Inow_Ipast_Gpast;
	unsigned long * F_Ipast_Jpast_Gpast;
	unsigned long * F_Inow_Ipast_Jpast_Gpast;
	gsl_vector_int * IndexMultipliers_Ipast_Gpast;
	gsl_vector_int * IndexMultipliers_Inow_Ipast_Gpast;
	gsl_vector_int * IndexMultipliers_Ipast_Jpast_Gpast;
	gsl_vector_int * IndexMultipliers_Inow_Ipast_Jpast_Gpast;
	unsigned long Tspace, Sspace;

	gsl_vector_int * vec_Full;
	gsl_vector_int * vec_Full_Bins;
	
	gsl_vector_int_view vec_Inow;
	gsl_vector_int_view vec_Ipast;
	gsl_vector_int_view vec_Jpast;
	// here the conditioning signal is fixed to order 1
	gsl_vector_int_view vec_Gpast;
	// int vec_G;
	
  rawdata **xdata;
	rawdata *xglobal;
#ifndef SEPARATED_OUTPUT
  double **xresult;
	long double Hxx, Hxxy;
#else
  double ***xresult;
	long double *Hxx, *Hxxy;
#endif
	double *xtest;
	rawdata *xtestD;

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
		sim.get("rawdatabins",rawdatabins);
    sim.get("AutoBinNumberQ",AutoBinNumberQ,false);
    if(!AutoBinNumberQ) sim.get("bins",bins);
		sim.get("globalbins",globalbins);
		sim.get("samples",samples);
		sim.get("StartSampleIndex",StartSampleIndex,1);
		sim.get("EndSampleIndex",EndSampleIndex,samples-1);
		sim.get("EqualSampleNumberQ",EqualSampleNumberQ,false);
		sim.get("MaxSampleNumberPerBin",MaxSampleNumberPerBin,-1);

		sim.get("noise",std_noise,-1.);
		sim.get("appliedscaling",input_scaling,1.);
		sim.get("cutoff",cutoff,-1.);
		sim.get("tauF",tauF);
		sim.get("OverrideRescalingQ",OverrideRescalingQ,false);
		sim.get("HighPassFilterQ",HighPassFilterQ,false);
		sim.get("InstantFeedbackTermQ",InstantFeedbackTermQ,false);
#ifdef ENABLE_ADAPTIVE_BINNING_AT_COMPILE_TIME
		sim.get("AdaptiveBinningQ",AdaptiveBinningQ,false);
		assert((!AdaptiveBinningQ)||(bins==2));
#endif
		sim.get("saturation",fluorescence_saturation,-1.);
		sim.get("IncludeGlobalSignalQ",IncludeGlobalSignalQ,true);
		assert(IncludeGlobalSignalQ ^ (globalbins==1));
		sim.get("GenerateGlobalFromFilteredDataQ",GenerateGlobalFromFilteredDataQ,false);
    sim.get("AutoConditioningLevelQ",AutoConditioningLevelQ,false);
		if(!AutoConditioningLevelQ)
		{
		  sim.get("GlobalConditioningLevel",GlobalConditioningLevel,-1.);
  		if (GlobalConditioningLevel>0) assert(globalbins==2); // for now, to keep it simple
  	}
		
		sim.get("SourceMarkovOrder",SourceMarkovOrder,1);
		assert(SourceMarkovOrder>0);
		// speed test:
		if (SourceMarkovOrder==1) SingleSourceMarkovOrder = true;
		else SingleSourceMarkovOrder = false;
		sim.get("TargetMarkovOrder",TargetMarkovOrder,1);
		assert(TargetMarkovOrder>0);
		// speed test:
		if (TargetMarkovOrder==1) SingleTargetMarkovOrder = true;
		else SingleTargetMarkovOrder = false;
		
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

		sim.get("ContinueOnErrorQ",ContinueOnErrorQ,false);

		// initialize random number generator
    // gsl_rng_env_setup();
    // GSLrandom = gsl_rng_alloc(gsl_rng_default);
    // gsl_rng_set(GSLrandom, 1234);
		
		AvailableSamples = NULL;
	};

	void execute(Sim& sim)
	{
	  sim.io <<"------ transferentropy-sim:extended ------ olav, Tue 12 Oct 2010 ------"<<Endl;
	  time_t start, middle, end;

	  time(&start);
	  sim.io <<"start: "<<ctime(&start)<<Endl;
	  sim.io <<"output file: "<<outputfile_results_name<<Endl;
		// Gespeichert wird später - hier nur Test, ob das Zielverzeichnis existiert
		write_parameters();
		skip_the_rest = false;

	  sim.io <<"allocating memory..."<<Endl;
		try {
      // xdata = NULL;
      // xdata = new rawdata*[size];
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
		
      // xglobal = new rawdata[samples];
#ifdef SEPARATED_OUTPUT
			// for testing
			Hxx = new long double[globalbins];
			Hxxy = new long double[globalbins];
#endif
			AvailableSamples = new unsigned long[globalbins];
      sim.io <<" -> done."<<Endl;			
			
      double** xdatadouble;
      if(inputfile_name!="")
        sim.io <<"input file: \""<<inputfile_name<<"\""<<Endl;
      else {
        sim.io <<"input files:"<<Endl;
        sim.io <<"- spike indices: \""<<spikeindexfile_name<<"\""<<Endl;
        sim.io <<"- spike times: \""<<spiketimesfile_name<<"\""<<Endl;
      }
      
			if(inputfile_name=="") {
        sim.io <<"loading data and generating time series from spike data..."<<Endl;
        xdatadouble = generate_time_series_from_spike_data(spiketimesfile_name, spikeindexfile_name, size, int(round(tauF)), samples, FluorescenceModel, std_noise, fluorescence_saturation, cutoff, DeltaCalciumOnAP, tauCa);
      }
      else {
        sim.io <<"loading data from binary file..."<<Endl;
        xdatadouble = load_time_series_from_binary_file(inputfile_name, size, samples, input_scaling, OverrideRescalingQ, std_noise, fluorescence_saturation, cutoff, sim);
      }
      sim.io <<" -> done."<<Endl;
      
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
      
      if(!GenerateGlobalFromFilteredDataQ) {
        sim.io <<"generating discretized global signal..."<<Endl;
        xglobal = generate_discretized_global_time_series(xdatadouble, size, samples, globalbins, GlobalConditioningLevel, AvailableSamples, StartSampleIndex, EndSampleIndex, sim);
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
			Tspace = Sspace = 1;
			for (int i=1; i<=TargetMarkovOrder; i++)
				Tspace *= bins;
			for (int i=1; i<=SourceMarkovOrder; i++)
				Sspace *= bins;
			F_Ipast_Gpast = NULL;
			F_Inow_Ipast_Gpast = NULL;
			F_Ipast_Jpast_Gpast = NULL;
			F_Inow_Ipast_Jpast_Gpast = NULL;
			F_Ipast_Gpast = new unsigned long[Tspace*globalbins];
			F_Inow_Ipast_Gpast = new unsigned long[bins*Tspace*globalbins];
			F_Ipast_Jpast_Gpast = new unsigned long[Tspace*Sspace*globalbins];
			F_Inow_Ipast_Jpast_Gpast = new unsigned long[bins*Tspace*Sspace*globalbins];
		
			// This is overall iterator that will be mapped onto array indices later:
			vec_Full = NULL;
			vec_Full_Bins = NULL;
			vec_Full = gsl_vector_int_alloc(1+TargetMarkovOrder+SourceMarkovOrder+1);
			gsl_vector_int_set_zero(vec_Full);
			vec_Full_Bins = gsl_vector_int_alloc(1+TargetMarkovOrder+SourceMarkovOrder+1);
			gsl_vector_int_set_all(vec_Full_Bins,bins);
			gsl_vector_int_set(vec_Full_Bins,1+TargetMarkovOrder+SourceMarkovOrder,globalbins);
		
			// Initialize views to have better access to the full iterator elements:
			vec_Inow = gsl_vector_int_subvector(vec_Full,0,1);
			vec_Ipast = gsl_vector_int_subvector(vec_Full,1,TargetMarkovOrder);
			vec_Jpast = gsl_vector_int_subvector(vec_Full,1+TargetMarkovOrder,SourceMarkovOrder);
			vec_Gpast = gsl_vector_int_subvector(vec_Full,1+TargetMarkovOrder+SourceMarkovOrder,1);
           
      if(GenerateGlobalFromFilteredDataQ) {
        sim.io <<"generating discretized global signal..."<<Endl;
        xglobal = generate_discretized_global_time_series(xdatadouble, size, samples, globalbins, GlobalConditioningLevel, AvailableSamples, StartSampleIndex, EndSampleIndex, sim);
        sim.io <<" -> done."<<Endl;
      }

      sim.io <<"discretizing local time series..."<<Endl;
      xdata = generate_discretized_version_of_time_series(xdatadouble, size, samples, bins);
      // and, since double version of time series is not used any more...
      free_time_series_memory(xdatadouble, size);
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
    // sim.io <<" -> done."<<Endl;
	
		// cout <<"testing discretization:"<<endl;
		// xtest = new double[100];
		// xtestD = new rawdata[100];
		// for(int i=0;i<100;i++)
		// 	xtest[i] = i;
		// discretize(xtest,xtestD,smallest(xtest,100),largest(xtest,100),100,bins);
		// for(int i=0;i<100;i++)
		// 	cout <<i<<": "<<xtest[i]<<" -> "<<(int)xtestD[i]<<endl;
		// exit(0);
		
		// cout <<"testing GSL vector class:"<<endl;
		// bool runningI = true;
		// while(runningI)
		// {
		// 	SimplePrintGSLVector(vec_Full);
		// 	runningI = OneStepAhead_FullIterator();
		// }
		
		if (!skip_the_rest) {
		SetUpMultidimArrayIndexMultipliers();
		
		// cout <<"--- TESTING ---"<<endl;
		// memset(F_Ipast_Gpast,0,sizeof(unsigned long)*Tspace*globalbins);
		// memset(F_Inow_Ipast_Gpast,0,sizeof(unsigned long)*bins*Tspace*globalbins);
		// memset(F_Ipast_Jpast_Gpast,0,sizeof(unsigned long)*Tspace*Sspace*globalbins);
		// memset(F_Inow_Ipast_Jpast_Gpast,0,sizeof(unsigned long)*bins*Tspace*Sspace*globalbins);
		// gsl_vector_int_set_all(vec_Full,1);
		// F_Ipast_Gpast[CounterArrayIndex(COUNTARRAY_IPAST_GPAST)]++;
		// F_Inow_Ipast_Gpast[CounterArrayIndex(COUNTARRAY_INOW_IPAST_GPAST)]++;
		// F_Ipast_Jpast_Gpast[CounterArrayIndex(COUNTARRAY_IPAST_JPAST_GPAST)]++;
		// F_Inow_Ipast_Jpast_Gpast[CounterArrayIndex(COUNTARRAY_INOW_IPAST_JPAST_GPAST)]++;
		// 
		// gsl_vector_int_set_zero(vec_Full);
		// bool runningIt = true;
		// unsigned long ig, iig, ijg,iijg;
		// while (runningIt)
		// {
		// 	// SimplePrintGSLVector(vec_Full);
		// 	ig = F_Ipast_Gpast[CounterArrayIndex(COUNTARRAY_IPAST_GPAST)];
		// 	iig = F_Inow_Ipast_Gpast[CounterArrayIndex(COUNTARRAY_INOW_IPAST_GPAST)];
		// 	ijg = F_Ipast_Jpast_Gpast[CounterArrayIndex(COUNTARRAY_IPAST_JPAST_GPAST)];
		// 	iijg = F_Inow_Ipast_Jpast_Gpast[CounterArrayIndex(COUNTARRAY_INOW_IPAST_JPAST_GPAST)];
		// 	// cout <<"CounterArrayIndex: "<<CounterArrayIndex(COUNTARRAY_INOW_IPAST_JPAST_GPAST)<<endl;
		// 	// if(ig>0){ cout <<"ig > 0 at: "; SimplePrintFullIterator(); cout <<" - CI:"<<CounterArrayIndex(COUNTARRAY_IPAST_GPAST)<<endl;}
		// 	// if(iig>0){ cout <<"iig > 0 at: "; SimplePrintFullIterator(); cout <<" - CI:"<<CounterArrayIndex(COUNTARRAY_INOW_IPAST_GPAST)<<endl;}
		// 	// if(ijg>0){ cout <<"ijg > 0 at: "; SimplePrintFullIterator(); cout <<" - CI:"<<CounterArrayIndex(COUNTARRAY_IPAST_JPAST_GPAST)<<endl;}
		// 	if(iijg>0){ cout <<"iijg > 0 at: "; SimplePrintFullIterator(); cout <<" - CI:"<<CounterArrayIndex(COUNTARRAY_INOW_IPAST_JPAST_GPAST)<<endl;}
		// 
		// 	runningIt = OneStepAhead_FullIterator();
		// }
		// 
		// exit(1);
		
#ifdef ENABLE_MODEL_FROM_SPIKES_AT_COMPILE_TIME
		generate_data_from_spikes();
#endif

    //    sim.io <<"loading data and adding noise (std "<<std_noise<<") and generating global signal... "<<Endl;
    //    load_data();
    // if (EqualSampleNumberQ) sim.io <<" (enforcing equal sample number per global bin)"<<Endl;
    //    sim.io <<" -> done."<<Endl;
	
	  // main loop:
		sim.io <<"set-up: "<<size<<" nodes, ";
		sim.io <<EndSampleIndex-StartSampleIndex+1<<" out of "<<samples<<" samples, ";
		sim.io <<bins<<" bins, "<<globalbins<<" globalbins"<<Endl;
		sim.io <<"set-up: Markov order of source/target/conditioning: "<<SourceMarkovOrder<<"/"<<TargetMarkovOrder<<"/1"<<Endl;
#ifdef SEPARATED_OUTPUT
		sim.io <<"set-up: separated output (globalbin)"<<Endl;
#endif
	
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
			if ((!status_already_displayed)&&((ii>=size/3)||(middle-start>30.)))
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
	      	xresult[jj][ii] = TransferEntropy(xdata[ii], xdata[jj]);
#else
	      	TransferEntropySeparated(xdata[ii], xdata[jj], ii, jj);
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
			write_multidim_result(xresult,globalbins);
#else
		  write_result(xresult);
#endif
			write_parameters();

			// free allocated memory
      // gsl_rng_free(GSLrandom);
			gsl_vector_int_free(vec_Full);
			gsl_vector_int_free(vec_Full_Bins);
    	gsl_vector_int_free(IndexMultipliers_Ipast_Gpast);
    	gsl_vector_int_free(IndexMultipliers_Inow_Ipast_Gpast);
    	gsl_vector_int_free(IndexMultipliers_Ipast_Jpast_Gpast);
    	gsl_vector_int_free(IndexMultipliers_Inow_Ipast_Jpast_Gpast);
		}

		try {
#ifdef SEPARATED_OUTPUT
		for (int x=0; x<size; x++)
		{
  		for (int x2=0; x2<size; x2++)
  			delete[] xresult[x][x2];
			delete[] xresult[x];
  	}
    delete[] Hxx;
    delete[] Hxxy;
#endif
		delete[] xresult;
		
		if (AvailableSamples != NULL) delete[] AvailableSamples;

		if (F_Ipast_Gpast != NULL) delete[] F_Ipast_Gpast;
		if (F_Inow_Ipast_Gpast != NULL) delete[] F_Inow_Ipast_Gpast;
		if (F_Ipast_Jpast_Gpast != NULL) delete[] F_Ipast_Jpast_Gpast;
		if (F_Inow_Ipast_Jpast_Gpast != NULL) delete[] F_Inow_Ipast_Jpast_Gpast;

    // for (int x=0; x<bins; x++)
    //  if (xdata[x] != NULL) delete[] xdata[x];
    // if (xdata != NULL) delete[] xdata;
    // if (xglobal != NULL) delete[] xglobal;
    free_time_series_memory(xdata,size);
    free_time_series_memory(xglobal);
		}
		catch(...) {};
		
		sim.io <<"End of Kernel (iteration="<<(sim.iteration())<<")"<<Endl;
	};
	

#ifdef SEPARATED_OUTPUT
	void TransferEntropySeparated(rawdata *arrayI, rawdata *arrayJ, int I, int J)
#else
	double TransferEntropy(rawdata *arrayI, rawdata *arrayJ)
#endif
	{
		// see for reference:
		//      Gourevitch und Eggermont. Evaluating Information Transfer Between Auditory
		//      Cortical Neurons. Journal of Neurophysiology (2007) vol. 97 (3) pp. 2533
		// We are looking at the information flow of array1 ("J") -> array2 ("I")
	
		// clear memory
#ifdef SEPARATED_OUTPUT
		memset(Hxx, 0, globalbins*sizeof(long double));
		memset(Hxxy, 0, globalbins*sizeof(long double));
#else
	  double result = 0.0;
		Hxx = Hxxy = 0.0;
#endif
		memset(F_Ipast_Gpast,0,sizeof(unsigned long)*Tspace*globalbins);
		memset(F_Inow_Ipast_Gpast,0,sizeof(unsigned long)*bins*Tspace*globalbins);
		memset(F_Ipast_Jpast_Gpast,0,sizeof(unsigned long)*Tspace*Sspace*globalbins);
		memset(F_Inow_Ipast_Jpast_Gpast,0,sizeof(unsigned long)*bins*Tspace*Sspace*globalbins);
	
	  // extract probabilities (actually number of occurrence)
		unsigned long const JShift = 0 + 1*InstantFeedbackTermQ;
		assert(StartSampleIndex >= max(TargetMarkovOrder,SourceMarkovOrder));
		for (unsigned long t=StartSampleIndex; t<EndSampleIndex; t++)
	  {
			// prepare the index vector vec_Full via the vector views
			gsl_vector_int_set(&vec_Inow.vector,0,arrayI[t]);
			// for (int i=0; i<TargetMarkovOrder; i++)
			// 	gsl_vector_int_set(&vec_Ipast.vector,i,arrayI[t-1+JShift-i]);
			// for (int i=0; i<SourceMarkovOrder; i++)
			// 	gsl_vector_int_set(&vec_Jpast.vector,i,arrayJ[t-1+JShift-i]);
			// speed test:
			if (SingleTargetMarkovOrder)
				gsl_vector_int_set(&vec_Ipast.vector,0,arrayI[t-1]);
			else for (int i=0; i<TargetMarkovOrder; i++)
				gsl_vector_int_set(&vec_Ipast.vector,i,arrayI[t-1-i]);
				
			if (SingleSourceMarkovOrder)
				gsl_vector_int_set(&vec_Jpast.vector,0,arrayJ[t-1+JShift]);
			else for (int i=0; i<SourceMarkovOrder; i++)
				gsl_vector_int_set(&vec_Jpast.vector,i,arrayJ[t-1+JShift-i]);
			gsl_vector_int_set(&vec_Gpast.vector,0,xglobal[t-1+JShift]);
			
			// add counts to arrays
			if (xglobal[t-1+JShift]<globalbins) { // for EqualSampleNumberQ flag
				F_Ipast_Gpast[CounterArrayIndex(COUNTARRAY_IPAST_GPAST)]++;
				F_Inow_Ipast_Gpast[CounterArrayIndex(COUNTARRAY_INOW_IPAST_GPAST)]++;
				F_Ipast_Jpast_Gpast[CounterArrayIndex(COUNTARRAY_IPAST_JPAST_GPAST)]++;
				F_Inow_Ipast_Jpast_Gpast[CounterArrayIndex(COUNTARRAY_INOW_IPAST_JPAST_GPAST)]++;
			}

			// DEBUG: test countings
			// if (t<50)
			// {
			// 	cout <<"t = "<<t<<", F_Full: ";
			// 	SimplePrintFullIterator(false);
			// 	cout <<", F_Inow_Ipast_Jpast_Gpast = "<<F_Inow_Ipast_Jpast_Gpast[CounterArrayIndex(COUNTARRAY_INOW_IPAST_JPAST_GPAST)]<<endl;
			// 
			// }
				
		}

		// DEBUG: test countings
		// gsl_vector_int_set_zero(vec_Full);
		// bool runningIt2 = true;
		// while (runningIt2)
		// {
		// 	SimplePrintFullIterator(false);
		// 	cout <<", F_Inow_Ipast_Jpast_Gpast = "<<F_Inow_Ipast_Jpast_Gpast[CounterArrayIndex(COUNTARRAY_INOW_IPAST_JPAST_GPAST)]<<endl;
		// 	
		// 	runningIt2 = OneStepAhead_FullIterator();
		// }
		// exit(0);

		// Here is some space for elaborate debiasing... :-)
		
		// Calculate transfer entropy from plug-in estimator:
		gsl_vector_int_set_zero(vec_Full);
		bool runningIt = true;
		unsigned long ig, iig, ijg,iijg;
		long double igd, iigd, ijgd,iijgd;
		while (runningIt)
		{
			// SimplePrintGSLVector(vec_Full);
			ig = F_Ipast_Gpast[CounterArrayIndex(COUNTARRAY_IPAST_GPAST)];
			igd = (long double)ig;
			iig = F_Inow_Ipast_Gpast[CounterArrayIndex(COUNTARRAY_INOW_IPAST_GPAST)];
			if (iig!=0)
			{
				iigd = (long double)iig;
				ijg = F_Ipast_Jpast_Gpast[CounterArrayIndex(COUNTARRAY_IPAST_JPAST_GPAST)];
				ijgd = (long double)ijg;
				iijg = F_Inow_Ipast_Jpast_Gpast[CounterArrayIndex(COUNTARRAY_INOW_IPAST_JPAST_GPAST)];
				iijgd = (long double)iijg;
			
				// a.) calculate Hxx:
				if (gsl_vector_int_isnull(&vec_Jpast.vector)) // to avoid counting the same term multiple times
				{
					if (iig!=0)
#ifdef SEPARATED_OUTPUT
						Hxx[gsl_vector_int_get(&vec_Gpast.vector,0)] -= iigd/AvailableSamples[gsl_vector_int_get(&vec_Gpast.vector,0)] * log(iigd/igd);
#else
						Hxx -= iigd/AvailableSamples[gsl_vector_int_get(&vec_Gpast.vector,0)] * log(iigd/igd);				
#endif
				}
			
				// b.) calculate Hxxy:
				if (iijg!=0)
#ifdef SEPARATED_OUTPUT
					Hxxy[gsl_vector_int_get(&vec_Gpast.vector,0)] -= iijgd/AvailableSamples[gsl_vector_int_get(&vec_Gpast.vector,0)] * log(iijgd/ijgd);				
#else
					Hxxy -= iijgd/AvailableSamples[gsl_vector_int_get(&vec_Gpast.vector,0)] * log(iijgd/ijgd);				
#endif
			}
			runningIt = OneStepAhead_FullIterator();
		}
		
#ifdef SEPARATED_OUTPUT
		for (rawdata g=0; g<globalbins; g++)
		{
			Hxx[g] /= log(2);
			Hxxy[g] /= log(2);
			xresult[J][I][g] = double(Hxx[g] - Hxxy[g]);
		}
#else
	  return double((Hxx - Hxxy)/log(2));
#endif

		// DEBUG
		// cout <<endl;
		// for (char g=0; g<globalbins; g++)
		// {
		// 	cout <<"Hxx="<<Hxx[g]<<", Hxxy="<<Hxxy[g]<<endl;
		// 	cout <<"-> result for g="<<int(g)<<": "<<xresult[J][I][g]<<endl;
		// }
		// exit(0);
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
		fileout1 <<"executable->teextendedsim";
		fileout1 <<", iteration->"<<iteration;
		
		fileout1 <<", size->"<<size;
		fileout1 <<", rawdatabins->"<<rawdatabins;
		fileout1 <<", bins->"<<bins;
		fileout1 <<", cutoff->"<<cutoff;
		fileout1 <<", globalbins->"<<globalbins;
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

		fileout1 <<", noise->"<<std_noise;
		fileout1 <<", tauF->"<<tauF;
		fileout1 <<", OverrideRescalingQ->"<<bool2textMX(OverrideRescalingQ);
		fileout1 <<", HighPassFilterQ->"<<bool2textMX(HighPassFilterQ);
		fileout1 <<", InstantFeedbackTermQ->"<<bool2textMX(InstantFeedbackTermQ);
#ifdef ENABLE_ADAPTIVE_BINNING_AT_COMPILE_TIME
		fileout1 <<", AdaptiveBinningQ->"<<bool2textMX(AdaptiveBinning)Q;
#endif
		fileout1 <<", saturation->"<<fluorescence_saturation;
		fileout1 <<", IncludeGlobalSignalQ->"<<bool2textMX(IncludeGlobalSignalQ);
		fileout1 <<", GenerateGlobalFromFilteredDataQ->"<<bool2textMX(GenerateGlobalFromFilteredDataQ);
		fileout1 <<", GlobalConditioningLevel->"<<GlobalConditioningLevel;
		fileout1 <<", TargetMarkovOrder->"<<TargetMarkovOrder;
		fileout1 <<", SourceMarkovOrder->"<<SourceMarkovOrder;
		
		fileout1 <<", AutoBinNumberQ->"<<AutoBinNumberQ;
    fileout1 <<", AutoConditioningLevelQ->"<<AutoConditioningLevelQ;
    
		
		fileout1 <<", inputfile->\""<<inputfile_name<<"\"";
		fileout1 <<", outputfile->\""<<outputfile_results_name<<"\"";
#ifdef ENABLE_MODEL_FROM_SPIKES_AT_COMPILE_TIME
		fileout1 <<", spikeindexfile->\""<<spikeindexfile_name<<"\"";
		fileout1 <<", spiketimesfile->\""<<spiketimesfile_name<<"\"";
		fileout1 <<", FluorescenceModel->\""<<FluorescenceModel<<"\"";
#endif
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
		for(int i=0; i<=1+TargetMarkovOrder+SourceMarkovOrder+1; i++)
		{
			if (i==1+TargetMarkovOrder+SourceMarkovOrder+1) return false; // if we have reached the "maximum" value
			gsl_vector_int_set(vec_Full,i,gsl_vector_int_get(vec_Full,i)+1);
			
			if (gsl_vector_int_get(vec_Full,i) >= gsl_vector_int_get(vec_Full_Bins,i))
			{
				gsl_vector_int_set(vec_Full,i,0);
				if((i==1+TargetMarkovOrder+SourceMarkovOrder)&&(GlobalConditioningLevel>0.0))
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
		cout <<"Inow: "<<gsl_vector_int_get(&vec_Inow.vector,0);
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
	inline unsigned long CounterArrayIndex(char ArrayCode)
	{
		// This could probably be optimized by using smarter IndexMultipliers, such that one can use
		// the GSL BLAS scalar product...
		
		unsigned long result = 0;
		
		switch (ArrayCode)
		{
			case COUNTARRAY_IPAST_GPAST: // dim = TargetMarkovOrder + 1
				if (SingleTargetMarkovOrder)
					result += gsl_vector_int_get(IndexMultipliers_Ipast_Gpast,0)*gsl_vector_int_get(&vec_Ipast.vector,0);
				else for (int i=0; i<TargetMarkovOrder; i++)
					result += gsl_vector_int_get(IndexMultipliers_Ipast_Gpast,i)*gsl_vector_int_get(&vec_Ipast.vector,i);
				result += gsl_vector_int_get(IndexMultipliers_Ipast_Gpast,TargetMarkovOrder)*gsl_vector_int_get(&vec_Gpast.vector,0);
				// assert (result<Tspace*globalbins);
				break;
			
			case COUNTARRAY_INOW_IPAST_GPAST: // dim. = 1 + TargetMarkovOrder + 1
				// result += gsl_vector_int_get(IndexMultipliers_Inow_Ipast_Gpast,0)*gsl_vector_int_get(&vec_Inow.vector,0);
				// for (int i=0; i<TargetMarkovOrder; i++)
				// 	result += gsl_vector_int_get(IndexMultipliers_Inow_Ipast_Gpast,i+1)*gsl_vector_int_get(&vec_Ipast.vector,i);
				// hopefully much faster:
				for (int i=0; i<1+TargetMarkovOrder; i++)
					result += gsl_vector_int_get(IndexMultipliers_Inow_Ipast_Gpast,i)*gsl_vector_int_get(vec_Full,i);				
				result += gsl_vector_int_get(IndexMultipliers_Inow_Ipast_Gpast,TargetMarkovOrder+1)*gsl_vector_int_get(&vec_Gpast.vector,0);
				// assert (result<bins*Tspace*globalbins);
				break;
			
			case COUNTARRAY_IPAST_JPAST_GPAST: // dim = TargetMarkovOrder + SourceMarkovOrder + 1
				// for (int i=0; i<TargetMarkovOrder; i++)
				// 	result += gsl_vector_int_get(IndexMultipliers_Ipast_Jpast_Gpast,i)*gsl_vector_int_get(&vec_Ipast.vector,i);
				// for (int i=0; i<SourceMarkovOrder; i++)
				// 	result += gsl_vector_int_get(IndexMultipliers_Ipast_Jpast_Gpast,i+TargetMarkovOrder)*gsl_vector_int_get(&vec_Jpast.vector,i);
				// hopefully much faster:
				for (int i=1; i<1+TargetMarkovOrder+SourceMarkovOrder+1; i++)
					result += gsl_vector_int_get(IndexMultipliers_Ipast_Jpast_Gpast,i-1)*gsl_vector_int_get(vec_Full,i);
				// result += gsl_vector_int_get(IndexMultipliers_Ipast_Jpast_Gpast,TargetMarkovOrder+SourceMarkovOrder)*gsl_vector_int_get(&vec_Gpast.vector,0);
				// assert (result<Tspace*Sspace*globalbins);
				break;
			
			case COUNTARRAY_INOW_IPAST_JPAST_GPAST: // dim = 1 + TargetMarkovOrder + SourceMarkovOrder + 1
				// result += gsl_vector_int_get(IndexMultipliers_Inow_Ipast_Jpast_Gpast,0)*gsl_vector_int_get(&vec_Inow.vector,0);
				// for (int i=0; i<TargetMarkovOrder; i++)
				// 	result += gsl_vector_int_get(IndexMultipliers_Inow_Ipast_Jpast_Gpast,i+1)*gsl_vector_int_get(&vec_Ipast.vector,i);
				// for (int i=0; i<SourceMarkovOrder; i++)
				// 	result += gsl_vector_int_get(IndexMultipliers_Inow_Ipast_Jpast_Gpast,i+TargetMarkovOrder+1)*gsl_vector_int_get(&vec_Jpast.vector,i);
				// result += gsl_vector_int_get(IndexMultipliers_Inow_Ipast_Jpast_Gpast,1+TargetMarkovOrder+SourceMarkovOrder)*gsl_vector_int_get(&vec_Gpast.vector,0);
				// assert (result<bins*Tspace*Sspace*globalbins);
				// hopefully much faster:
				for (int i=0; i<1+TargetMarkovOrder+SourceMarkovOrder+1; i++)
					result += gsl_vector_int_get(IndexMultipliers_Inow_Ipast_Jpast_Gpast,i)*gsl_vector_int_get(vec_Full,i);
				break;
			
			default: cout <<endl<<"GetCounterArray: error, invalid ArrayCode."<<endl; exit(1);
		}
		
		// if (result>=bins*Tspace*Sspace*globalbins)
		// {
		// 	cout <<endl<<"GetCounterArray: error, invalid result (ArrayCode="<<int(ArrayCode)<<", result="<<result;
		// 	cout <<", vec_Full: ";
		// 	SimplePrintFullIterator();
		// 	cout <<")"<<endl;
		// 	exit(1);
		// }
		return result;
	}
	
	// These multipliers are calculated for optimization: the actual indices are then just
	// a linear summation of terms like (pre-calculated multiplier * current index)
	void SetUpMultidimArrayIndexMultipliers()
	{
		// here we assume equal binning for source and target terms
		
		// BinsPerDim is actually a slight repetion of the logic for vec_Full_Bins - maybe there is a more elegant way...
		
		// because (SourceMarkovOrder+TargetMarkovOrder+2) is the maximal length:
		gsl_vector_int* BinsPerDim = gsl_vector_int_alloc(SourceMarkovOrder+TargetMarkovOrder+2);
		
		// ------------------ IndexMultipliers_Ipast_Gpast:
		// cout <<"Setting up IndexMultipliers_Ipast_Gpast:"<<endl;
		IndexMultipliers_Ipast_Gpast = gsl_vector_int_alloc(TargetMarkovOrder+1);
		gsl_vector_int_set_all(IndexMultipliers_Ipast_Gpast,1);
		gsl_vector_int_set_all(BinsPerDim,1);
		
		for (int i=0; i<TargetMarkovOrder; i++)
			gsl_vector_int_set(BinsPerDim,i,bins);
		gsl_vector_int_set(BinsPerDim,TargetMarkovOrder,globalbins);
		
		for (int i=1; i<IndexMultipliers_Ipast_Gpast->size; i++)
			gsl_vector_int_set(IndexMultipliers_Ipast_Gpast, i, gsl_vector_int_get(BinsPerDim,i-1)*gsl_vector_int_get(IndexMultipliers_Ipast_Gpast, i-1));	
		
		// cout <<"BinsPerDim: ";
		// SimplePrintGSLVector(BinsPerDim);
		// cout <<"IndexMultipliers_Ipast_Gpast: ";
		// SimplePrintGSLVector(IndexMultipliers_Ipast_Gpast);

		// ------------------ IndexMultipliers_Inow_Ipast_Gpast:
		// cout <<"DEBUG: Setting up IndexMultipliers_Inow_Ipast_Gpast:"<<endl;
		IndexMultipliers_Inow_Ipast_Gpast = gsl_vector_int_alloc(TargetMarkovOrder+2);
		gsl_vector_int_set_all(IndexMultipliers_Inow_Ipast_Gpast,1);
		gsl_vector_int_set_all(BinsPerDim,1);
		
		for (int i=0; i<1+TargetMarkovOrder; i++)
			gsl_vector_int_set(BinsPerDim,i,bins);
		gsl_vector_int_set(BinsPerDim,TargetMarkovOrder+1,globalbins);
		
		for (int i=1; i<IndexMultipliers_Inow_Ipast_Gpast->size; i++)
			gsl_vector_int_set(IndexMultipliers_Inow_Ipast_Gpast, i, gsl_vector_int_get(BinsPerDim,i-1)*gsl_vector_int_get(IndexMultipliers_Inow_Ipast_Gpast, i-1));	
		
		// cout <<"BinsPerDim: ";
		// SimplePrintGSLVector(BinsPerDim);
		// cout <<"IndexMultipliers_Inow_Ipast_Gpast: ";
		// SimplePrintGSLVector(IndexMultipliers_Inow_Ipast_Gpast);

		// ------------------ IndexMultipliers_Ipast_Jpast_Gpast:
		// cout <<"Setting up IndexMultipliers_Ipast_Jpast_Gpast:"<<endl;
		IndexMultipliers_Ipast_Jpast_Gpast = gsl_vector_int_alloc(TargetMarkovOrder+SourceMarkovOrder+1);
		gsl_vector_int_set_all(IndexMultipliers_Ipast_Jpast_Gpast,1);
		gsl_vector_int_set_all(BinsPerDim,1);
		
		for (int i=0; i<TargetMarkovOrder+SourceMarkovOrder; i++)
			gsl_vector_int_set(BinsPerDim,i,bins);
		gsl_vector_int_set(BinsPerDim,TargetMarkovOrder+SourceMarkovOrder,globalbins);
		
		for (int i=1; i<IndexMultipliers_Ipast_Jpast_Gpast->size; i++)
			gsl_vector_int_set(IndexMultipliers_Ipast_Jpast_Gpast, i, gsl_vector_int_get(BinsPerDim,i-1)*gsl_vector_int_get(IndexMultipliers_Ipast_Jpast_Gpast, i-1));	
		
		// cout <<"BinsPerDim: ";
		// SimplePrintGSLVector(BinsPerDim);
		// cout <<"IndexMultipliers_Ipast_Jpast_Gpast: ";
		// SimplePrintGSLVector(IndexMultipliers_Ipast_Jpast_Gpast);
		
		// ------------------ IndexMultipliers_Inow_Ipast_Jpast_Gpast:
		// cout <<"DEBUG: Setting up IndexMultipliers_Inow_Ipast_Jpast_Gpast:"<<endl;
		IndexMultipliers_Inow_Ipast_Jpast_Gpast = gsl_vector_int_alloc(1+TargetMarkovOrder+SourceMarkovOrder+1);
		gsl_vector_int_set_all(IndexMultipliers_Inow_Ipast_Jpast_Gpast,1);
		gsl_vector_int_set_all(BinsPerDim,1);
		
		for (int i=0; i<1+TargetMarkovOrder+SourceMarkovOrder; i++)
			gsl_vector_int_set(BinsPerDim,i,bins);
		gsl_vector_int_set(BinsPerDim,1+TargetMarkovOrder+SourceMarkovOrder,globalbins);
		
		for (int i=1; i<IndexMultipliers_Inow_Ipast_Jpast_Gpast->size; i++)
			gsl_vector_int_set(IndexMultipliers_Inow_Ipast_Jpast_Gpast, i, gsl_vector_int_get(BinsPerDim,i-1)*gsl_vector_int_get(IndexMultipliers_Inow_Ipast_Jpast_Gpast, i-1));	
		
		// cout <<"- BinsPerDim: ";
		// SimplePrintGSLVector(BinsPerDim);
		// cout <<"- IndexMultipliers_Inow_Ipast_Jpast_Gpast: ";
		// SimplePrintGSLVector(IndexMultipliers_Inow_Ipast_Jpast_Gpast);
		gsl_vector_int_free(BinsPerDim);
	};
	
	unsigned long SimpleScalarProduct(gsl_vector_int* v1, gsl_vector_int* v2)
	{
		assert(v1->size == v2->size);
		unsigned long sum = 0;
		for(int i=0; i<v1->size; i++)
			sum += gsl_vector_int_get(v1,i) * gsl_vector_int_get(v2,i);
		return sum;
	}
};
