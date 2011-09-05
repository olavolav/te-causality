// calculate the transfer entropy between a numer of time series
// this is the extension to arbitrary Markov order of the source and target term
// including the binless estimators based on FLANN nearest neighbor search
// created by olav, So  3 Jul 2011 14:31:37 CEST

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
#include <gsl/gsl_sf_psi.h>

#include "../../Simulationen/olav.h"
#include "../../../Sonstiges/SimKernel/sim_main.h"
#include "../te-datainit.h"

#include <flann/flann.hpp>
#define FLANN_NEIGHBORS 10
#define FLANN_NEIGHBORS_IN_MARGINALIZED_SPACES 10000*FLANN_NEIGHBORS

#undef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE

// #define GSL_RANDOM_NUMBER_GENERATOR gsl_rng_default
#define GSL_RANDOM_NUMBER_GENERATOR gsl_rng_ranlxs2

#define CALCUATE_TRANSER_ENTROPY_ONLY_FOR_LOWEST_GLOBAL_BIN

#define REPORTS 25
#undef SHOW_DETAILED_PROGRESS

#define OUTPUTNUMBER_PRECISION 15
#define SEPARATED_OUTPUT

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

	bool IncludeGlobalSignalQ;
	bool GenerateGlobalFromFilteredDataQ;
	double GlobalConditioningLevel;
	int SourceMarkovOrder, TargetMarkovOrder;
	
	bool ContinueOnErrorQ;
	bool skip_the_rest;
	
  // bool AutoBinNumberQ;
  bool AutoConditioningLevelQ;
  
  gsl_rng* GSLrandom;

  // container vectors for binless estimators
#ifndef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
  flann::Matrix<double>** datasetFLANN_Ipast;
  // flann::Matrix<int>* indicesFLANN_Ipast; flann::Matrix<double>* distancesFLANN_Ipast;
  flann::Matrix<double>** datasetFLANN_Inow_Ipast;
  // flann::Matrix<int>* indicesFLANN_Inow_Ipast; flann::Matrix<double>* distancesFLANN_Inow_Ipast;
#endif
  flann::Matrix<double>** datasetFLANN_Inow;
  // flann::Matrix<int>* indicesFLANN_Inow; flann::Matrix<double>* distancesFLANN_Inow;
  flann::Matrix<double>** datasetFLANN_Ipast_Jpast;
  // flann::Matrix<int>* indicesFLANN_Ipast_Jpast; flann::Matrix<double>* distancesFLANN_Ipast_Jpast;
  flann::Matrix<double>** datasetFLANN_Inow_Ipast_Jpast;
  // flann::Matrix<int>* indicesFLANN_Inow_Ipast_Jpast; flann::Matrix<double>* distancesFLANN_Inow_Ipast_Jpast;
  
  flann::Matrix<int>* indicesFLANN_radiusSearch;
  flann::Matrix<double>* distancesFLANN_radiusSearch;
#ifndef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
  flann::Matrix<int>* indicesFLANN_Inow_Ipast;
  flann::Matrix<double>* distancesFLANN_Inow_Ipast;
#endif
  flann::Matrix<int>* indicesFLANN_Inow_Ipast_Jpast;
  flann::Matrix<double>* distancesFLANN_Inow_Ipast_Jpast;

#ifdef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
  long* shuffle_permutation;
#endif

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
		if(std_noise<0) {
      sim.io <<"Noise hast to be positive for the NN estimator to work. Aborting."<<Endl;
      exit(1);
		}
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
		if(!AutoConditioningLevelQ) {
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
    if(!((inputfile_name!="") ^ ((spikeindexfile_name!="")&&(spiketimesfile_name!="")))) {
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
		
#ifndef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
    datasetFLANN_Ipast = NULL;
    // indicesFLANN_Ipast = NULL; distancesFLANN_Ipast = NULL;
    datasetFLANN_Inow_Ipast = NULL;
    // indicesFLANN_Inow_Ipast = NULL; distancesFLANN_Inow_Ipast = NULL;
#endif
    datasetFLANN_Inow = NULL;
    // indicesFLANN_Inow = NULL; distancesFLANN_Inow = NULL;
    datasetFLANN_Ipast_Jpast = NULL;
    // indicesFLANN_Ipast_Jpast = NULL; distancesFLANN_Ipast_Jpast = NULL;
    datasetFLANN_Inow_Ipast_Jpast = NULL;
    // indicesFLANN_Inow_Ipast_Jpast = NULL; distancesFLANN_Inow_Ipast_Jpast = NULL;

    indicesFLANN_radiusSearch = NULL;
    distancesFLANN_radiusSearch = NULL;
#ifndef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
    indicesFLANN_Inow_Ipast = NULL;
    distancesFLANN_Inow_Ipast = NULL;
#endif
    indicesFLANN_Inow_Ipast_Jpast = NULL;
    distancesFLANN_Inow_Ipast_Jpast = NULL;

#ifdef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
    shuffle_permutation = NULL;
#endif

		AvailableSamples = NULL;
    xdatadouble = NULL;
    xglobal = NULL;
    xresult = NULL;
	};

	void execute(Sim& sim)
	{
	  sim.io <<"------ transferentropy-sim:binless ------ olav, Wed 08 Jun 2011 ------"<<Endl;
    // time_t start, middle, end;

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
		  for(int i=0; i<size; i++) {
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
		
			AvailableSamples = new unsigned long[globalbins];
						
		  // hack of medium ugliness to make it work without global signal
		  if(globalbins<=1) {
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
        xdatadouble = generate_time_series_from_spike_data(spiketimesfile_name, spikeindexfile_name, size, int(round(tauF)), samples, FluorescenceModel, std_noise, fluorescence_saturation, cutoff, DeltaCalciumOnAP, tauCa, GSLrandom, sim);
      }
      else {
        sim.io <<"loading data from binary file..."<<Endl;
        xdatadouble = load_time_series_from_binary_file(inputfile_name, size, samples, input_scaling, OverrideRescalingQ, std_noise, fluorescence_saturation, cutoff, GSLrandom, sim);
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

      // reserve memory for FLANN structures
      int dims;
#ifndef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
      datasetFLANN_Ipast = new flann::Matrix<double>*[globalbins];
      datasetFLANN_Inow_Ipast = new flann::Matrix<double>*[globalbins];
#endif
      datasetFLANN_Inow = new flann::Matrix<double>*[globalbins];
      datasetFLANN_Ipast_Jpast = new flann::Matrix<double>*[globalbins];
      datasetFLANN_Inow_Ipast_Jpast = new flann::Matrix<double>*[globalbins];

      for(rawdata g=0; g<globalbins; g++) {
#ifndef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
        dims = TargetMarkovOrder; // marginal space -> radius search
        datasetFLANN_Ipast[g] = \
          new flann::Matrix<double>(new double[AvailableSamples[g]*dims],AvailableSamples[g],dims);
          
        dims = 1+TargetMarkovOrder; // joint space -> k-th NN search
        datasetFLANN_Inow_Ipast[g] = \
          new flann::Matrix<double>(new double[AvailableSamples[g]*dims],AvailableSamples[g],dims);
        indicesFLANN_Inow_Ipast = new flann::Matrix<int>(new int[FLANN_NEIGHBORS_IN_MARGINALIZED_SPACES+1],1,FLANN_NEIGHBORS_IN_MARGINALIZED_SPACES+1);
#endif
        dims = 1; // marginal space -> radius search
        datasetFLANN_Inow[g] = \
          new flann::Matrix<double>(new double[AvailableSamples[g]*dims],AvailableSamples[g],dims);

        dims = TargetMarkovOrder+SourceMarkovOrder; // marginal space -> radius search
        datasetFLANN_Ipast_Jpast[g] = \
          new flann::Matrix<double>(new double[AvailableSamples[g]*dims],AvailableSamples[g],dims);
          
        dims = 1+TargetMarkovOrder+SourceMarkovOrder; // joint space -> k-th NN search
        datasetFLANN_Inow_Ipast_Jpast[g] = \
          new flann::Matrix<double>(new double[AvailableSamples[g]*dims],AvailableSamples[g],dims);
      }
      
      // indices and distances for searches in marginalized spaces (radius search with only 1 query point)
      indicesFLANN_radiusSearch = new flann::Matrix<int>(new int[FLANN_NEIGHBORS_IN_MARGINALIZED_SPACES+1],1,FLANN_NEIGHBORS_IN_MARGINALIZED_SPACES+1);
      distancesFLANN_radiusSearch = new flann::Matrix<double>(new double[FLANN_NEIGHBORS_IN_MARGINALIZED_SPACES+1],1,FLANN_NEIGHBORS_IN_MARGINALIZED_SPACES+1);
      // indices and distances for searches in joint spaces (kNN search with many query points)
#ifndef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
      indicesFLANN_Inow_Ipast = new flann::Matrix<int>(new int[effectivesamples*(FLANN_NEIGHBORS+1)],effectivesamples,FLANN_NEIGHBORS+1);
      distancesFLANN_Inow_Ipast = new flann::Matrix<double>(new double[effectivesamples*(FLANN_NEIGHBORS+1)],effectivesamples,FLANN_NEIGHBORS+1);
#endif
      indicesFLANN_Inow_Ipast_Jpast = new flann::Matrix<int>(new int[effectivesamples*(FLANN_NEIGHBORS+1)],effectivesamples,FLANN_NEIGHBORS+1);
      distancesFLANN_Inow_Ipast_Jpast = new flann::Matrix<double>(new double[effectivesamples*(FLANN_NEIGHBORS+1)],effectivesamples,FLANN_NEIGHBORS+1);

#ifdef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
      shuffle_permutation = generate_random_geometric_permutation(samples,globalbins,xglobal,5,GSLrandom);
#endif
		}
		catch(...) {
			sim.io <<"Error: could not reserve enough memory!"<<Endl;
			if(!ContinueOnErrorQ) exit(1);
			else {
				sim.io <<"Error handling: ContinueOnErrorQ flag set, proceeding..."<<Endl;
				skip_the_rest = true;
			}
		}
    // sim.io <<" -> done."<<Endl;
			
		if (!skip_the_rest) {		
  	  // ----------------------------------------- main loop: start -----------------------------------------
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
  			if ((!status_already_displayed)&&((ii>=size/3)||((middle-start>30.)&&(ii>0)))) { 
  				sim.io <<" (after "<<ii<<" nodes: elapsed "<<sec2string(difftime(middle,start)) \
  					<<", ETA "<<ETAstring(ii,size,difftime(middle,start))<<")"<<Endl;
  				status_already_displayed = true;
  			}
#endif
  	    for(int jj=0; jj<size; jj++)
  	    {
          // cout <<"debug: launching TE for "<<jj<<" -> "<<ii<<" ..."<<endl;
  	      if (ii != jj) {
#ifndef SEPARATED_OUTPUT
  	      	xresult[jj][ii] = TransferEntropyKraskov(xdatadouble[ii], xdatadouble[jj]);
#else
  	      	TransferEntropyKraskov(xdatadouble[ii], xdatadouble[jj], ii, jj);
#endif
  	      }
  	      // else xresult[ii][jj] = 0.0;
  	    }
  	  }
#ifndef SHOW_DETAILED_PROGRESS
	  sim.io <<" -> done."<<Endl;
#endif
    }
  	// ----------------------------------------- main loop: end -----------------------------------------
	  time(&end);
	  sim.io <<"end: "<<ctime(&end)<<Endl;
	  sim.io <<"runtime: "<<sec2string(difftime(end,start))<<Endl;

		// cout <<"TE terms: "<<terms_sum<<", of those zero: "<<terms_zero<<" ("<<int(double(terms_zero)*100/terms_sum)<<"%)"<<endl;
	};
	
	void finalize(Sim& sim) {
		if(!skip_the_rest) {
#ifdef SEPARATED_OUTPUT
			write_multidim_result(xresult,globalbins);
#else
		  write_result(xresult);
#endif
			write_parameters();

			// free allocated memory
      gsl_rng_free(GSLrandom);  
#ifdef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
      if (shuffle_permutation != NULL ) delete[] shuffle_permutation;
#endif
      
      for(rawdata g=0; g<globalbins; g++) {
#ifndef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
        delete[] datasetFLANN_Ipast[g]->data;
        delete[] datasetFLANN_Inow_Ipast[g]->data;
#endif
        delete[] datasetFLANN_Ipast_Jpast[g]->data;
        delete[] datasetFLANN_Inow_Ipast_Jpast[g]->data;
      }
#ifndef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
      if (datasetFLANN_Ipast != NULL) delete[] datasetFLANN_Ipast;
      if (datasetFLANN_Inow_Ipast != NULL) delete[] datasetFLANN_Inow_Ipast;
#endif
      if (datasetFLANN_Ipast_Jpast != NULL) delete[] datasetFLANN_Ipast_Jpast;
      if (datasetFLANN_Inow_Ipast_Jpast != NULL) delete[] datasetFLANN_Inow_Ipast_Jpast;

      delete[] indicesFLANN_radiusSearch->data;
      delete[] distancesFLANN_radiusSearch->data;
#ifndef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
      delete[] indicesFLANN_Inow_Ipast->data;
      delete[] distancesFLANN_Inow_Ipast->data;
#endif
      delete[] indicesFLANN_Inow_Ipast_Jpast->data;
      delete[] distancesFLANN_Inow_Ipast_Jpast->data;

  	}

		try {
#ifdef SEPARATED_OUTPUT
  		for (int x=0; x<size; x++) {
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


	std::string bool2textMX(bool value) {
		if (value) return "True";
		else return "False";
	};

	void write_parameters() {
		char* name = new char[outputfile_pars_name.length()+1];
		strcpy(name,outputfile_pars_name.c_str());
		ofstream fileout1(name);
		delete[] name;
		if (fileout1 == NULL) {
	  	cerr <<endl<<"error: cannot open output file!"<<endl;
	  	exit(1);
	  }

		fileout1.precision(6);
		fileout1 <<"{";
		fileout1 <<"executable->tebinlesssim";
#ifdef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
		fileout1 <<", NormalizationViaShuffling->True";
#else
		fileout1 <<", NormalizationViaShuffling->False";
#endif
		fileout1 <<", iteration->"<<iteration;
		time(&end);
    fileout1 <<", ExecutionTime->"<<sec2string(difftime(end,start));
		
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
		for (int i=0; i<globalbins; i++) {
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

	void write_result(double **array) {
		char* name = new char[outputfile_results_name.length()+1];
		strcpy(name,outputfile_results_name.c_str());
		ofstream fileout1(name);
		delete[] name;
		if (fileout1 == NULL) {
	  	cerr <<endl<<"error: cannot open output file!"<<endl;
	  	exit(1);
	  }	  

		fileout1.precision(OUTPUTNUMBER_PRECISION);
		fileout1 <<fixed;
	  fileout1 <<"{";
	  for(int j=0; j<size; j++) {
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

	void write_multidim_result(double ***array, unsigned int dimens) {
		char* name = new char[outputfile_results_name.length()+1];
		strcpy(name,outputfile_results_name.c_str());
		ofstream fileout1(name);
		delete[] name;
		if (fileout1 == NULL) {
	  	cerr <<endl<<"error: cannot open output file!"<<endl;
	  	exit(1);
	  }	  

		fileout1.precision(OUTPUTNUMBER_PRECISION);
		fileout1 <<fixed;
	  fileout1 <<"{";
	  for(unsigned long j=0; j<size; j++) {
	  	if(j>0) fileout1<<",";
	  	fileout1 <<"{";
	    for(unsigned long i=0; i<size; i++) {
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

	  cout <<"Transfer entropy matrix saved."<<endl;
	};
  
  long double MututalInformationKraskovFLANN(flann::Matrix<double>* datasetA, flann::Matrix<double>* datasetB, flann::Matrix<double>* datasetJoint, flann::Matrix<int>* indicesJoint, flann::Matrix<double>* distancesJoint, const long samples_here) {
    // reference:
    // Kraskov, A., Stögbauer, H., & Grassberger, P. (2004). Estimating mutual information. Physical Review E
    // (here we implement the first algorithm, see equation 8 in the paper)
    long double result = 0.;
    long neighborsA = 0;
    long neighborsB = 0;
    double epsilon;
    
    // print some dataset elements
    // for(int row=0; row<datasetA->rows; row++) {
    //   for(int col=0; col<datasetA->cols; col++) {
    //     cout <<"debug: d["<<row<<"]["<<col<<"] = "<<(*datasetA)[row][col]<<endl;
    //     // cout <<"debug: set dataset["<<row<<"]["<<col<<"] to "<<(*dataset)[row][col]<<endl;
    //   }
    // }
    
    // 0.) Set up hierarchical k-means clustering for NN search later
    // (actually, we only need the maximum norm for the joint space, other norm might be faster...)
    flann::Index<flann::MaxDistance<double> > indexA(*datasetA, flann::KMeansIndexParams());
    // flann::Index<flann::MaxDistance<double> > indexA(*datasetA, flann::AutotunedIndexParams());
    indexA.buildIndex();
    // exit(1);
    
    flann::Index<flann::MaxDistance<double> > indexB(*datasetB, flann::KMeansIndexParams());
    // default parameters from test file: ... and then: 32, 11, flann::FLANN_CENTERS_RANDOM, 0.2));
    indexB.buildIndex();
    flann::Index<flann::MaxDistance<double> > indexJoint(*datasetJoint, flann::KMeansIndexParams());
    indexJoint.buildIndex();
    
    // 1.) Search for distance eps to k-th neighbor to s in joint space
    // cout <<"debug: k-search ..."<<endl;
    indexJoint.knnSearch(*datasetJoint,*indicesJoint,*distancesJoint,FLANN_NEIGHBORS+1,flann::SearchParams(-1));
    
    // 2.) Search for the number of neighbors in the spaces A and B within a "radius" eps
    flann::Matrix<double>* qPointA = new flann::Matrix<double>(new double[1*datasetA->cols],1,datasetA->cols);
    flann::Matrix<double>* qPointB = new flann::Matrix<double>(new double[1*datasetB->cols],1,datasetB->cols);
    for(long s=0; s<samples_here; s++) {
      epsilon = (*distancesJoint)[s][FLANN_NEIGHBORS]; // here no +1 of course 
      // cout <<"debug: radius-search #1 ..."<<endl;
      // cout <<"bla: "<<(qPointA->cols)<<", "<<indexA.veclen()<<endl;
      for(long i=0; i<datasetA->cols; i++) {
        (*qPointA)[0][i] = (*datasetA)[s][i];
      }
      neighborsA = indexA.radiusSearch(*qPointA,*indicesFLANN_radiusSearch,*distancesFLANN_radiusSearch,epsilon/2., flann::SearchParams(-1));
      if(neighborsA == FLANN_NEIGHBORS_IN_MARGINALIZED_SPACES+1) {
        cout <<"Warning in MututalInformationKraskovFLANN: number of neighbors (radius, #A) may exceed memory limit: "<<endl;
      }
      // cout <<"debug: radius-search #2 ..."<<endl;
      for(long i=0; i<datasetB->cols; i++) {
        (*qPointB)[0][i] = (*datasetB)[s][i];
      }
      neighborsB = indexB.radiusSearch(*qPointB,*indicesFLANN_radiusSearch,*distancesFLANN_radiusSearch,epsilon/2.,flann::SearchParams(-1));
      if(neighborsB == FLANN_NEIGHBORS_IN_MARGINALIZED_SPACES+1) {
        cout <<"Warning in MututalInformationKraskovFLANN: number of neighbors (radius, #B) may exceed memory limit"<<endl;
      }
      // cout <<"debug: epsilon = "<<epsilon<<", neighborsA = "<<neighborsA<<", neighborsB = "<<neighborsB<<endl;
      result += gsl_sf_psi(neighborsA+1) + gsl_sf_psi(neighborsB+1);
    }
    
    // 3.) Compute the Mutual Information (add first and last term)
    result = gsl_sf_psi(FLANN_NEIGHBORS) + result/samples_here + gsl_sf_psi(samples_here);
    
    // cout <<"debug: MututalInformationKraskovFLANN = "<<result<<endl;
    return result;
  };


#ifdef SEPARATED_OUTPUT
	void TransferEntropyKraskov(double *arrayI, double *arrayJ, int I, int J)
#else
	double TransferEntropyKraskov(double *arrayI, double *arrayJ)
#endif
  {
#ifndef SEPARATED_OUTPUT
    xresult = 0.;
#endif
		const long JShift = (const long)InstantFeedbackTermQ;
    long t_sample;
    long double result = 0.;
    
    for(rawdata g=0; g<globalbins; g++) {
      // cout <<"debug: TransferEntropyKraskov launching globalbin "<<int(g)<<" ..."<<endl;
      if(AvailableSamples[g] > 1) {
        // 1.) set input structures
        t_sample = 0;
        for(long t=StartSampleIndex; t<=EndSampleIndex; t++) {
          if(xglobal[t] == g) {
            // set Inow
#ifndef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
            (*(datasetFLANN_Inow_Ipast[g]))[t_sample][0] = arrayI[t];
#endif
            (*(datasetFLANN_Inow[g]))[t_sample][0] = arrayI[t];
            (*(datasetFLANN_Inow_Ipast_Jpast[g]))[t_sample][0] = arrayI[t];
            // set Ipast
            for (int k=0; k<TargetMarkovOrder; k++) {
#ifndef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE
              (*(datasetFLANN_Ipast[g]))[t_sample][k] = arrayI[t-(k+1)];
              (*(datasetFLANN_Inow_Ipast[g]))[t_sample][1+k] = arrayI[t-(k+1)];
#endif
              (*(datasetFLANN_Ipast_Jpast[g]))[t_sample][k] = arrayI[t-(k+1)];
              (*(datasetFLANN_Inow_Ipast_Jpast[g]))[t_sample][1+k] = arrayI[t-(k+1)];
            }
            // set Jpast
            for (int l=0; l<SourceMarkovOrder; l++) {
              (*(datasetFLANN_Ipast_Jpast[g]))[t_sample][TargetMarkovOrder+l] = arrayJ[t-(l+1)+JShift];
              (*(datasetFLANN_Inow_Ipast_Jpast[g]))[t_sample][1+TargetMarkovOrder+l] = arrayJ[t-(l+1)+JShift];
            }
            t_sample++;
          }
        }
        assert(AvailableSamples[g] == t_sample);
        
        // test datasetFLANN_Ipast_Jpast:
        // for(long tt=0; tt<AvailableSamples[g]; tt++) {
        //   for (int l=0; l<TargetMarkovOrder+SourceMarkovOrder; l++) {
        //     assert((*(datasetFLANN_Ipast_Jpast[g]))[tt][l]<100.);
        // }}

        // 2.) calculate MI terms (based on correctly ordered source data)
        result += 
          MututalInformationKraskovFLANN(datasetFLANN_Ipast_Jpast[g], datasetFLANN_Inow[g], datasetFLANN_Inow_Ipast_Jpast[g], indicesFLANN_Inow_Ipast_Jpast, distancesFLANN_Inow_Ipast_Jpast, AvailableSamples[g]);
        // cout <<"debug in TransferEntropyKraskov: result1 =  "<<int(g)<<" ..."<<endl;
        
#ifndef NORMALIZE_TRANSFER_ENTROPY_ESTIMATE    
        result -= 
          MututalInformationKraskovFLANN(datasetFLANN_Ipast[g], datasetFLANN_Inow[g], datasetFLANN_Inow_Ipast[g], indicesFLANN_Inow_Ipast, distancesFLANN_Inow_Ipast, AvailableSamples[g]);

#else

        // 3.) shuffle source data
        t_sample = 0;
        for(long t=StartSampleIndex; t<=EndSampleIndex; t++) {
          if(xglobal[t] == g) {
            // set Jpast (because this is the only part that we shuffle)
            for (int l=0; l<SourceMarkovOrder; l++) {
              (*datasetFLANN_Ipast_Jpast[g])[t_sample][TargetMarkovOrder+l] = arrayJ[shuffle_permutation[t-(l+1)+JShift]];
              (*datasetFLANN_Inow_Ipast_Jpast[g])[t_sample][1+TargetMarkovOrder+l] = arrayJ[shuffle_permutation[t-(l+1)+JShift]];
            }
            t_sample++;
          }
        }
      
        // 4.) calculate differential entropy terms (based on shuffeled source data)
        result -= 
          MututalInformationKraskovFLANN(datasetFLANN_Ipast_Jpast[g], datasetFLANN_Inow[g], datasetFLANN_Inow_Ipast_Jpast[g], indicesFLANN_Inow_Ipast_Jpast, distancesFLANN_Inow_Ipast_Jpast, AvailableSamples[g]);
#endif

        // cout <<"debug: TE = "<<result<<endl;

#ifdef SEPARATED_OUTPUT
  	    xresult[J][I][g] = double(result);
#endif
      }
      else {
#ifdef SEPARATED_OUTPUT
  	    xresult[J][I][g] = 0.;
#endif
      }
#ifdef CALCUATE_TRANSER_ENTROPY_ONLY_FOR_LOWEST_GLOBAL_BIN
      break;
#endif
    } // end globalbins loop
  
#ifndef SEPARATED_OUTPUT
    return result;
#endif
  };

};
