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
#include <gsl/gsl_vector_double.h>

#include "../../Simulationen/olav.h"
#include "../../../Sonstiges/SimKernel/sim_main.h"
#include "../te-datainit.h"

// #ifndef INLINE
// #define INLINE extern inline
// #endif

#define REPORTS 25
// #define SHOW_DETAILED_PROGRESS

#define OUTPUTNUMBER_PRECISION 15
#define SEPARATED_OUTPUT

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
  // unsigned int bins
  unsigned int globalbins;
	// unsigned long mag der SimKernel irgendwie nicht?
	long samples, effectivesamples;
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
	
  // bool AutoBinNumberQ;
  bool AutoConditioningLevelQ;

  // container vectors for binless estimators
  // gsl_vector** Sample_Ipast;
  // gsl_vector** Sample_Inow_Ipast;
  gsl_vector** Sample_Ipast_Jpast;
  gsl_vector** Sample_Inow_Ipast_Jpast;
  gsl_rng* GSLshuffler;
  long* shuffle_permutation;
	
  // rawdata **xdata;
  double **xdatadouble;
	rawdata *xglobal;
#ifndef SEPARATED_OUTPUT
  double **xresult;
#else
  double ***xresult;
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
    // bins = 0;
    // sim.get("AutoBinNumberQ",AutoBinNumberQ,false);
    // if(!AutoBinNumberQ) sim.get("bins",bins);
		sim.get("globalbins",globalbins,1);
    if(globalbins<1) globalbins=1;
		sim.get("samples",samples);
		sim.get("StartSampleIndex",StartSampleIndex,1);
    assert(StartSampleIndex>=1 && StartSampleIndex<samples);
		sim.get("EndSampleIndex",EndSampleIndex,samples-1);
    assert(EndSampleIndex>=1 && EndSampleIndex<samples);
    effectivesamples = EndSampleIndex-StartSampleIndex+1;
		sim.get("EqualSampleNumberQ",EqualSampleNumberQ,false);
		sim.get("MaxSampleNumberPerBin",MaxSampleNumberPerBin,-1);

		sim.get("noise",std_noise,-1.);
		sim.get("appliedscaling",input_scaling,1.);
		sim.get("cutoff",cutoff,-1.);
		sim.get("tauF",tauF);
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
    GSLshuffler = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(GSLshuffler, 1234);
		
    Sample_Ipast_Jpast = NULL;
    Sample_Inow_Ipast_Jpast = NULL;
    shuffle_permutation = NULL;

		AvailableSamples = NULL;
    xdatadouble = NULL;
    xglobal = NULL;
    xresult = NULL;
	};

	void execute(Sim& sim)
	{
	  sim.io <<"------ transferentropy-sim:binless ------ olav, Wed 08 Jun 2011 ------"<<Endl;
	  time_t start, middle, end;

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
		
// #ifdef SEPARATED_OUTPUT
//      // for testing
//      Hxx = new long double[globalbins];
//      Hxxy = new long double[globalbins];
// #endif
			AvailableSamples = new unsigned long[globalbins];
			
      Sample_Ipast_Jpast = new gsl_vector*[effectivesamples];
      Sample_Inow_Ipast_Jpast = new gsl_vector*[effectivesamples];
      for(long t=0; t<effectivesamples; t++) {
        Sample_Ipast_Jpast[t] = gsl_vector_alloc(TargetMarkovOrder+SourceMarkovOrder);
        Sample_Inow_Ipast_Jpast[t] = gsl_vector_alloc(1+TargetMarkovOrder+SourceMarkovOrder);
      }
      shuffle_permutation = new long[samples]; // here we potentially need the whole range of samples
			
		  // hack of medium ugliness to make it work without global signal
		  if(globalbins<=1)
		  {
        // sim.io <<"debug: xglobal hack."<<Endl;
        xglobal = new rawdata[samples];
        memset(xglobal, 0, samples*sizeof(rawdata));
        AvailableSamples[0] = EndSampleIndex-StartSampleIndex+1;
		  }
      sim.io <<" -> done."<<Endl;			
			
      // double** xdatadouble = NULL;
      if(inputfile_name!="")
        sim.io <<"input file: \""<<inputfile_name<<"\""<<Endl;
      else {
        sim.io <<"input files:"<<Endl;
        sim.io <<"- spike indices: \""<<spikeindexfile_name<<"\""<<Endl;
        sim.io <<"- spike times: \""<<spiketimesfile_name<<"\""<<Endl;
      }
      
			if(inputfile_name=="") {
        sim.io <<"loading data and generating time series from spike data..."<<Endl;
        xdatadouble = generate_time_series_from_spike_data(spiketimesfile_name, spikeindexfile_name, size, int(round(tauF)), samples, FluorescenceModel, std_noise, fluorescence_saturation, cutoff, DeltaCalciumOnAP, tauCa, sim);
      }
      else {
        sim.io <<"loading data from binary file..."<<Endl;
        xdatadouble = load_time_series_from_binary_file(inputfile_name, size, samples, input_scaling, OverrideRescalingQ, std_noise, fluorescence_saturation, cutoff, sim);
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
        xglobal = generate_discretized_global_time_series(xdatadouble, size, samples, globalbins, GlobalConditioningLevel, AvailableSamples, StartSampleIndex, EndSampleIndex, sim);
        sim.io <<" -> done."<<Endl;
      }
      
      if(HighPassFilterQ) {
        sim.io <<"applying high-pass filter to time series..."<<Endl;
        apply_high_pass_filter_to_time_series(xdatadouble, size, samples);
        sim.io <<" -> done."<<Endl;
      }
      
      if((globalbins>1)&&(GenerateGlobalFromFilteredDataQ)) {
        sim.io <<"generating discretized global signal..."<<Endl;
        xglobal = generate_discretized_global_time_series(xdatadouble, size, samples, globalbins, GlobalConditioningLevel, AvailableSamples, StartSampleIndex, EndSampleIndex, sim);
        sim.io <<" -> done."<<Endl;
      }

      // sim.io <<"discretizing local time series..."<<Endl;
      // xdata = generate_discretized_version_of_time_series(xdatadouble, size, samples, bins);
      // and, since double version of time series is not used any more...
      // if(xdatadouble!=NULL) free_time_series_memory(xdatadouble, size);
      // sim.io <<" -> done."<<Endl;

      // cout <<"DEBUG: subset of discretized global signal: ";
      // display_subset(xglobal);
      
      // generate random permutation (once for all)
      for(long t=0; t<samples; t++)
        shuffle_permutation[t] = t;
      gsl_ran_shuffle(GSLshuffler,&shuffle_permutation[0],samples,sizeof(long));
      
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
			
		if (!skip_the_rest) {		
	
	  // main loop:
		sim.io <<"set-up: "<<size<<" nodes, ";
		sim.io <<EndSampleIndex-StartSampleIndex+1<<" out of "<<samples<<" samples, ";
    // sim.io <<bins<<" bins, "
		sim.io <<globalbins<<" globalbins"<<Endl;
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
	      	xresult[jj][ii] = DifferentialTransferEntropy(xdatadouble[ii], xdatadouble[jj]);
#else
	      	DifferentialTransferEntropy(xdatadouble[ii], xdatadouble[jj], ii, jj);
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
      gsl_rng_free(GSLshuffler);
      
      for(long t=0; t<effectivesamples; t++) {
        gsl_vector_free(Sample_Ipast_Jpast[t]);
        gsl_vector_free(Sample_Inow_Ipast_Jpast[t]);
      }
      if (Sample_Ipast_Jpast != NULL) delete[] Sample_Ipast_Jpast;
      if (Sample_Inow_Ipast_Jpast != NULL) delete[] Sample_Inow_Ipast_Jpast;
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
		
		if (AvailableSamples != NULL) delete[] AvailableSamples;

    if(xdatadouble != NULL) free_time_series_memory(xdatadouble,size);
    if(xglobal != NULL) free_time_series_memory(xglobal);
		}
		catch(...) {};
		
		sim.io <<"End of Kernel (iteration="<<(sim.iteration())<<")"<<Endl;
	};

	
#ifdef SEPARATED_OUTPUT
	void DifferentialTransferEntropy(double *arrayI, double *arrayJ, int I, int J)
#else
	double DifferentialTransferEntropy(double *arrayI, double *arrayJ)
#endif
  {
		long const JShift = (long const)InstantFeedbackTermQ;
    long t_sample;
    long double H_Ipast_Jpast, H_Inow_Ipast_Jpast;
    
    for(rawdata g=0; g<globalbins; g++)
    {    
      // 1.) set GSL-Vectors as input for the differential entropy calculation
      t_sample = 0;
      for(long t=StartSampleIndex; t<=EndSampleIndex; t++) {
        if(xglobal[t] == g) {
          // set Inow
          gsl_vector_set(Sample_Inow_Ipast_Jpast[t_sample],0,arrayI[t]);
          // set Ipast
          for (int k=0; k<TargetMarkovOrder; k++) {
            gsl_vector_set(Sample_Ipast_Jpast[t_sample],k,arrayI[t-(k+1)]);
            gsl_vector_set(Sample_Inow_Ipast_Jpast[t_sample],1+k,arrayI[t-(k+1)]);
          }
          // set Jpast
          for (int l=0; l<SourceMarkovOrder; l++) {
            gsl_vector_set(Sample_Ipast_Jpast[t_sample],TargetMarkovOrder+l,arrayJ[t-(l+1)+JShift]);
            gsl_vector_set(Sample_Inow_Ipast_Jpast[t_sample],1+TargetMarkovOrder+l,arrayJ[t-(l+1)-JShift]);
          }
          t_sample++;
        }
      }
      // cout <<"DEBUG: AvailableSamples["<<int(g)<<"] = "<<AvailableSamples[g];
      // cout <<", t_sample = "<<t_sample<<endl;
      assert(AvailableSamples[g] == t_sample);
    
      // 2.) calculate differential entropy terms (based on correctly ordered source data)
      H_Ipast_Jpast = DifferentialEntropy(Sample_Ipast_Jpast,TargetMarkovOrder+SourceMarkovOrder,t_sample);
      cout <<"DEBUG_1: H_Ipast_Jpast = "<<H_Ipast_Jpast<<endl;
      H_Inow_Ipast_Jpast = DifferentialEntropy(Sample_Inow_Ipast_Jpast,1+TargetMarkovOrder+SourceMarkovOrder,t_sample);
      cout <<"DEBUG_1: H_Inow_Ipast_Jpast = "<<H_Inow_Ipast_Jpast<<endl;
    
      // 3a.) generate random permutation (skipped, we use the same for all)
      // long shuffle_permutation[samples]; // here we potentially need the whole range of samples
      // for(long t=0; t<samples; t++)
      //   shuffle_permutation[t] = t;
      // gsl_ran_shuffle(GSLshuffler,shuffle_permutation,samples,sizeof(long));
      
      // 3b.) shuffle source data
      for(long t=StartSampleIndex; t<=EndSampleIndex; t++) {
        t_sample = t - StartSampleIndex;
        // set Jpast (because this is the only part that we shuffle)
        for (int l=0; l<SourceMarkovOrder; l++) {
          gsl_vector_set(Sample_Ipast_Jpast[t_sample], TargetMarkovOrder+l, arrayJ[shuffle_permutation[t-(l+1)+JShift]]);
          gsl_vector_set(Sample_Inow_Ipast_Jpast[t_sample], 1+TargetMarkovOrder+l, arrayJ[shuffle_permutation[t-(l+1)-JShift]]);
        }
      }
      
      // 4.) calculate differential entropy terms (based on shuffeled source data)
      H_Ipast_Jpast -= DifferentialEntropy(Sample_Ipast_Jpast,TargetMarkovOrder+SourceMarkovOrder,t_sample);
      cout <<"DEBUG_2: H_Ipast_Jpast = "<<H_Ipast_Jpast<<endl;
      H_Inow_Ipast_Jpast -= DifferentialEntropy(Sample_Inow_Ipast_Jpast,1+TargetMarkovOrder+SourceMarkovOrder,t_sample);
      cout <<"DEBUG_2: H_Inow_Ipast_Jpast = "<<H_Inow_Ipast_Jpast<<endl;
    
      // 5.) calculate result

#ifdef SEPARATED_OUTPUT
  	  xresult[J][I][g] = double(H_Inow_Ipast_Jpast-H_Ipast_Jpast);
#else
      xresult += double(H_Inow_Ipast_Jpast-H_Ipast_Jpast);
#endif
    } // end globalbins loop
  
#ifndef SEPARATED_OUTPUT
    return xresult;
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
		fileout1 <<"executable->teextendedsim";
		fileout1 <<", iteration->"<<iteration;
		
		fileout1 <<", size->"<<size;
    // fileout1 <<", bins->"<<bins;
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
};
