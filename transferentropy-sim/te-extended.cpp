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

// #ifndef INLINE
// #define INLINE extern inline
// #endif

#define REPORTS 25
#undef SHOW_DETAILED_PROGRESS

#define OUTPUTNUMBER_PRECISION 15
#define SEPARATED_OUTPUT

#undef ENABLE_MODEL_FROM_SPIKES_AT_COMPILE_TIME
#undef ENABLE_ADAPTIVE_BINNING_AT_COMPILE_TIME
#define FMODEL_SPIKECOUNT 1
#define FMODEL_LEOGANG_RAWCALCIUM 2
#define FMODEL_LEOGANG_INTEGRATED_RAWCALCIUM 3
#define FMODEL_ERROR -1

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
	unsigned int word_length;
	double std_noise;
	string inputfile_name;
	string outputfile_results_name;
	string outputfile_pars_name;
#ifdef ENABLE_MODEL_FROM_SPIKES_AT_COMPILE_TIME
	string spikeindexfile_name, spiketimesfile_name;
	string FluorescenceModel;
	int intFluorescenceModel;
#endif
	gsl_rng* GSLrandom;
	double input_scaling;
	double cutoff;
	double tauF;
	double fluorescence_saturation;
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
		sim.get("bins",bins);
		sim.get("globalbins",globalbins);
		sim.get("samples",samples);
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
		sim.get("GenerateGlobalFromFilteredDataQ",GenerateGlobalFromFilteredDataQ,false);
		sim.get("GlobalConditioningLevel",GlobalConditioningLevel,-1.);
		if (GlobalConditioningLevel>0) assert(globalbins==2); // for now, to keep it simple
		
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
#ifdef ENABLE_MODEL_FROM_SPIKES_AT_COMPILE_TIME
		sim.get("spikeindexfile",spikeindexfile_name,"");
		sim.get("spiketimesfile",spiketimesfile_name,"");
		sim.get("FluorescenceModel",FluorescenceModel,"");
		intFluorescenceModel = FMODEL_ERROR;
		if (FluorescenceModel=="SpikeCount")
			intFluorescenceModel = FMODEL_SPIKECOUNT;
		if (FluorescenceModel=="Leogang-rawCalcium")
			intFluorescenceModel = FMODEL_LEOGANG_RAWCALCIUM;
		if (FluorescenceModel=="Leogang-Integrated-rawCalcium")
			intFluorescenceModel = FMODEL_LEOGANG_INTEGRATED_RAWCALCIUM;
		// assert(intFluorescenceModel != FMODEL_ERROR);
#endif

		// initialize random number generator
		gsl_rng_env_setup();
		GSLrandom = gsl_rng_alloc(gsl_rng_default);
		gsl_rng_set(GSLrandom, 1234);
	};

	void execute(Sim& sim)
	{
	  sim.io <<"------ transferentropy-sim:extended ------ olav, Tue 12 Oct 2010 ------"<<Endl;
	  time_t start, end;
#ifdef SHOW_DETAILED_PROGRESS
	  time_t middle;
#endif
	  unsigned long totaltrials, completedtrials;

	  time(&start);
	  sim.io <<"start: "<<ctime(&start)<<Endl;
	  // sim.io <<"running on host: "<<flush;
	  // system("hostname");
	  // sim.io <<"current directory: "<<flush;
	  // system("pwd");

	  sim.io <<"input file: "<<inputfile_name<<Endl;
	  sim.io <<"output file: "<<outputfile_results_name<<Endl;

	  sim.io <<"allocating memory..."<<Endl;
	  xdata = new rawdata*[size];
#ifndef SEPARATED_OUTPUT
	  xresult = new double*[size];
#else
	  xresult = new double**[size];
#endif
	  for(int i=0; i<size; i++)
	  {
	    xdata[i] = new rawdata[samples];
	    memset(xdata[i], 0, samples*sizeof(rawdata));
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
	
		Tspace = Sspace = 1;
		for (int i=1; i<=TargetMarkovOrder; i++)
			Tspace *= bins;
		for (int i=1; i<=SourceMarkovOrder; i++)
			Sspace *= bins;
		F_Ipast_Gpast = new unsigned long[Tspace*globalbins];
		F_Inow_Ipast_Gpast = new unsigned long[bins*Tspace*globalbins];
		F_Ipast_Jpast_Gpast = new unsigned long[Tspace*Sspace*globalbins];
		F_Inow_Ipast_Jpast_Gpast = new unsigned long[bins*Tspace*Sspace*globalbins];
		
		// This is overall iterator that can and will be mapped onto array indices later:
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
		
		xglobal = new rawdata[samples];
#ifdef SEPARATED_OUTPUT
		// for testing
		Hxx = new long double[globalbins];
		Hxxy = new long double[globalbins];
#endif
	  sim.io <<" done."<<Endl;
	
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

	  sim.io <<"loading data and adding noise (std "<<std_noise<<") and generating global signal... "<<Endl;
	  load_data();
	  // generate_global();
	  sim.io <<" done."<<Endl;
	
	  // main loop:
	  totaltrials = size*(size-1);
	  sim.io <<"set-up: "<<size<<" neurons, "<<samples<<" samples, "<<bins<<" bins, "<<globalbins<<" globalbins"<<Endl;
		sim.io <<"set-up: Markov order of source/target/conditioning: "<<SourceMarkovOrder<<"/"<<TargetMarkovOrder<<"/1"<<Endl;
#ifdef SEPARATED_OUTPUT
		sim.io <<"set-up: separated output (globalbin)"<<Endl;
#endif
	  completedtrials = 0;
		// unsigned long long terms_sum = 0;
		// unsigned long long terms_zero = 0;
	
#ifdef SHOW_DETAILED_PROGRESS
			sim.io <<"running...";
#else
	  	sim.io <<"running "<<Endl;
#endif

	  for(int ii=0; ii<size; ii++)
	  {
#ifdef SHOW_DETAILED_PROGRESS
	  	status(ii,REPORTS,size);
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
					completedtrials++;
	      }
	      // else xresult[ii][jj] = 0.0;
	    }
	  }
#ifndef SHOW_DETAILED_PROGRESS
	  cout <<" done."<<endl;
#endif

	  time(&end);
	  sim.io <<"end: "<<ctime(&end)<<Endl;
	  sim.io <<"runtime: "<<sec2string(difftime(end,start))<<Endl;

		// cout <<"TE terms: "<<terms_sum<<", of those zero: "<<terms_zero<<" ("<<int(double(terms_zero)*100/terms_sum)<<"%)"<<endl; 
	};
	
	void finalize(Sim& sim)
	{
#ifdef SEPARATED_OUTPUT
		write_multidim_result(xresult,globalbins);
#else
	  write_result(xresult);
#endif
		save_parameters();

		// free allocated memory
		gsl_rng_free(GSLrandom);
		gsl_vector_int_free(vec_Full);
		gsl_vector_int_free(vec_Full_Bins);

#ifdef SEPARATED_OUTPUT
		for (int x=0; x<globalbins; x++)
			delete[] xresult[x];
#endif
		delete[] xresult;

		delete[] F_Ipast_Gpast;
		delete[] F_Inow_Ipast_Gpast;
		delete[] F_Ipast_Jpast_Gpast;
		delete[] F_Inow_Ipast_Jpast_Gpast;

		delete[] xdata;
		delete[] xglobal;
				
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
		unsigned long const word_length = max(TargetMarkovOrder, SourceMarkovOrder);
		for (unsigned long t=word_length; t<samples; t++)
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
			F_Ipast_Gpast[CounterArrayIndex(COUNTARRAY_IPAST_GPAST)]++;
			F_Inow_Ipast_Gpast[CounterArrayIndex(COUNTARRAY_INOW_IPAST_GPAST)]++;
			F_Ipast_Jpast_Gpast[CounterArrayIndex(COUNTARRAY_IPAST_JPAST_GPAST)]++;
			F_Inow_Ipast_Jpast_Gpast[CounterArrayIndex(COUNTARRAY_INOW_IPAST_JPAST_GPAST)]++;

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
						Hxx[gsl_vector_int_get(&vec_Gpast.vector,0)] -= iigd/(samples-word_length) * log(iigd/igd);
#else
						Hxx -= iigd/(samples-word_length) * log(iigd/igd);				
#endif
				}
			
				// b.) calculate Hxxy:
				if (iijg!=0)
#ifdef SEPARATED_OUTPUT
					Hxxy[gsl_vector_int_get(&vec_Gpast.vector,0)] -= iijgd/(samples-word_length) * log(iijgd/ijgd);				
#else
					Hxxy -= iijgd/(samples-word_length) * log(iijgd/ijgd);				
#endif
			}
			runningIt = OneStepAhead_FullIterator();
		}
		
#ifdef SEPARATED_OUTPUT
		for (char g=0; g<globalbins; g++)
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

	
	void load_data()
	{
		char* name = new char[inputfile_name.length()+1];
		strcpy(name,inputfile_name.c_str());
	  ifstream binaryfile(name, ios::binary);
		delete[] name;
		
		char* temparray = new char[samples];
		double* tempdoublearray = new double[samples];
		memset(tempdoublearray, 0, samples*sizeof(double));
		// double xtemp;
		memset(xglobal, 0, samples*sizeof(rawdata));
		double* xglobaltemp = new double[samples];
		memset(xglobaltemp, 0, samples*sizeof(double));

	  if (binaryfile == NULL)
	  {
	  	cout <<endl<<"error: cannot find input file!"<<endl;
	  	exit(1);
	  }
	
		// test file length
		binaryfile.seekg(0,ios::end);
		if(long(binaryfile.tellg()) != size*samples)
		{
	  	cout <<endl<<"error: file length of input does not match given parameters!"<<endl;
	  	exit(1);
		}
		binaryfile.seekg(0,ios::beg);

		double* tempdoublearraycopy = new double[samples];
		
	  for(int j=0; j<size; j++)
	  {
	    binaryfile.read(temparray, samples);
	
			// OverrideRescalingQ = true
			// Dies ignoriert also "appliedscaling", "noise", "HighPassFilterQ" und "cutoff"
			// Therefore, "bins" takes the role of an upper cutoff
			if (OverrideRescalingQ)
		    for(long k=0; k<samples; k++)
					xdata[j][k] = temparray[k];					
			else     // OverrideRescalingQ = false
			{
		    for (long k=0; k<samples; k++)
				{
					// transform to unsigned notation
					tempdoublearray[k] = double(temparray[k]);
					if (temparray[k]<0) tempdoublearray[k] += 256.;
					
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
				
				// add to what later becomes the global signal - depending on the position of this block
				// relative to the HighPass block, the global signal is generated from the filtered or
				// unfiltered signal
				if(GenerateGlobalFromFilteredDataQ == false)
			    for(unsigned long k=0; k<samples; k++) xglobaltemp[k] += tempdoublearray[k];

				if(HighPassFilterQ) {
					// of course, this is just a difference signal, so not really filtered
					memcpy(tempdoublearraycopy,tempdoublearray,samples*sizeof(double));
					tempdoublearray[0] = 0.0;
			    for(long k=1; k<samples; k++)
						tempdoublearray[k] = tempdoublearraycopy[k] - tempdoublearraycopy[k-1];
				}

				if(GenerateGlobalFromFilteredDataQ == true)
			    for(long k=0; k<samples; k++) xglobaltemp[k] += tempdoublearray[k];
			
#ifndef ENABLE_ADAPTIVE_BINNING_AT_COMPILE_TIME
				discretize(tempdoublearray,xdata[j],bins);
#else
				if(!AdaptiveBinningQ) discretize(tempdoublearray,xdata[j],bins);
				else discretize2accordingtoStd(tempdoublearray,xdata[j]);
#endif
			}
	  }
	
		// cout <<endl;
		// for(int j=0; j<400; j++)
		// 	cout <<int(xdata[2][j])<<",";
		// cout <<endl;
		// exit(1);

		// generate global signal (rescaled to globalbins binning)
		if (IncludeGlobalSignalQ)
		{
			for (unsigned long t=0; t<samples; t++)
				xglobaltemp[t] /= size;
				
			if (GlobalConditioningLevel > 0.)
			{
				unsigned long below = 0;
				for (unsigned long t=0; t<samples; t++)
				{
					if (xglobaltemp[t] > GlobalConditioningLevel) xglobal[t] = 1;
					else
					{
						xglobal[t] = 0;
						below++;
					}
				}
				cout <<"global conditioning with level "<<GlobalConditioningLevel<<": "<<(100.*below)/samples<<"% are below threshold... "<<endl;
			}
			else discretize(xglobaltemp,xglobal,globalbins);
		}
		else
			for (unsigned long t=0; t<samples; t++) xglobaltemp[t] = 0;
	
		delete[] temparray;
		delete[] xglobaltemp;
		delete[] tempdoublearray;
	};
	
	void discretize(double* in, rawdata* out, unsigned int nr_bins)
	{
		discretize(in,out,smallest(in,samples),largest(in,samples),samples,nr_bins);
	};
	void discretize(double* in, rawdata* out, double min, double max, unsigned int nr_samples, unsigned int nr_bins)
	{
		// correct discretization according to 'te-test.nb'
		double xstepsize = (max-min)/nr_bins;
		// cout <<"max = "<<max<<endl;
		// cout <<"min = "<<min<<endl;
		// cout <<"bins here = "<<nr_bins<<endl;
		// cout <<"stepsize = "<<xstepsize<<endl;

		int xint;
		for (unsigned long t=0; t<nr_samples; t++)
		{
			assert(in[t]<=max);
			if (in[t]>=max) xint = nr_bins-1;
			else
			{
				if (in[t]<=min) xint = 0;
				else xint = (int)((in[t]-min)/xstepsize); // if(!(xint<nr_bins)) cout<<"!"<<xint<<","<<in[t]<<endl; }
			}
			if (xint >= nr_bins) xint = nr_bins-1; // need to have this for some silly numerical reason...
			
			assert(xint>=0);

			out[t] = (rawdata)xint;
			assert(out[t]<nr_bins);
		}
	};

#ifdef ENABLE_ADAPTIVE_BINNING_AT_COMPILE_TIME	
	void discretize2accordingtoStd(double* in, rawdata* out)
	{
		double splitheight = sqrt(gsl_stats_variance(in,1,samples));

		int xint;
		for (unsigned long t=0; t<samples; t++)
		{
			if (in[t]>splitheight) out[t] = 1;
			else out[t] = 0;
		}
	};
#endif
	
#ifdef ENABLE_MODEL_FROM_SPIKES_AT_COMPILE_TIME
	// void generate_data_from_spikes (unsigned int length, double** outputarray)
  void generate_data_from_spikes ()
	{
		cout <<"generating time series form spike data..."<<endl;
		// open files
		char* nameI = new char[spikeindexfile_name.length()+1];
		strcpy(nameI,inputfile_name.c_str());
	  ifstream binaryfileI(nameI, ios::binary);
		delete[] nameI;
		char* nameT = new char[spiketimesfile_name.length()+1];
		strcpy(nameT,inputfile_name.c_str());
	  ifstream binaryfileT(nameT, ios::binary);
		delete[] nameT;
		
		// determine file length
		binaryfileI.seekg(0,ios::end);
		const long nr_spikes = binaryfileI.tellg();
		cout <<"number of spikes in index file: "<<nr_spikes<<endl;
		binaryfileI.seekg(0,ios::beg);
		int* xindex = new int[nr_spikes];
		double* xtimes = new double[nr_spikes];
		
		// load spike data
		binaryfileI.read((char*)xindex, nr_spikes*sizeof(int));
		binaryfileT.read(reinterpret_cast<char*>(xtimes), nr_spikes*sizeof(double));
		
		// generate fluorescence data
		const int int_tauF = (int)round(tauF); // in ms
		unsigned long startindex = 1;
		unsigned long endindex = 0;
		unsigned long ttExactMS = 0;
		for (unsigned long ttExactMS=0; ttExactMS<100+0*samples; ttExactMS+=tauF)
		{
			// cout <<"xindex = "<<(long)xindex[t]<<", xtimes = "<<xtimes[t]<<endl;
			// assuming that the data import works...
			
			// determine starting and ending spike index of current frame
			while ((endindex+1<nr_spikes)&&(xtimes[endindex+1]<=ttExactMS+int_tauF))
				endindex++;
				
			for (int ii=0; ii<size; ii++)
			{
				switch (intFluorescenceModel)
				{
					case FMODEL_SPIKECOUNT:
						
						break;
					case FMODEL_LEOGANG_RAWCALCIUM:
						break;
					case FMODEL_LEOGANG_INTEGRATED_RAWCALCIUM:
						break;
					default:
						cout <<"error in generate_data_from_spikes: invalid fluorescence model";
						exit(1);
				}
			}
			
		}
		
		delete[] xindex;
		delete[] xtimes;
	}
	
	unsigned long count(int* array, unsigned long starti, unsigned long endi, int what)
	{
		unsigned long occur = 0;
		for (unsigned long i=starti; i<=endi; i++)
			if (array[i] == what) occur++;
		return occur;
	}
#endif
	
	double smallest(double* array, unsigned int length)
	{
		double min = array[0];
		for (unsigned int i=1; i<length; i++)
			if(array[i]<min) min = array[i];

		return min;
	};
	double largest(double* array, unsigned int length)
	{
		double max = array[0];
		for (unsigned int i=1; i<length; i++)
			if(array[i]>max) max = array[i];

		return max;
	};

	void save_parameters()
	{
		char* name = new char[outputfile_pars_name.length()+1];
		strcpy(name,outputfile_pars_name.c_str());
		ofstream fileout1(name);
		delete[] name;
		if (fileout1 == NULL)
	  {
	  	cerr <<endl<<"error: cannot open parameters output file!"<<endl;
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
		fileout1 <<", noise->"<<std_noise;
		fileout1 <<", tauF->"<<tauF;
		fileout1 <<", OverrideRescalingQ->"<<OverrideRescalingQ;
		fileout1 <<", HighPassFilterQ->"<<HighPassFilterQ;
		fileout1 <<", InstantFeedbackTermQ->"<<InstantFeedbackTermQ;
#ifdef ENABLE_ADAPTIVE_BINNING_AT_COMPILE_TIME
		fileout1 <<", AdaptiveBinningQ->"<<AdaptiveBinningQ;
#endif
		fileout1 <<", saturation->"<<fluorescence_saturation;
		fileout1 <<", IncludeGlobalSignalQ->"<<IncludeGlobalSignalQ;
		fileout1 <<", GenerateGlobalFromFilteredDataQ->"<<GenerateGlobalFromFilteredDataQ;
		fileout1 <<", GlobalConditioningLevel->"<<GlobalConditioningLevel;
		fileout1 <<", TargetMarkovOrder->"<<TargetMarkovOrder;
		fileout1 <<", SourceMarkovOrder->"<<SourceMarkovOrder;
		
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
	
	void SimplePrintGSLVector(gsl_vector_int* vec)
	{
		SimplePrintGSLVector(vec,true);
	};
	void SimplePrintGSLVector(gsl_vector_int* vec, bool newline)
	{
		for(int i=0; i<vec->size; i++)
			cout <<gsl_vector_int_get(vec,i)<<" ";
		if (newline) cout <<endl;
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
