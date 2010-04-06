// calculate the mutual information between a numer of time series
// created by olav, Mo 5 Apr 2010 00:34:21 CEST

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <ctime>
#include <cstring>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>

#include "../../Simulationen/olav.h"
#include "../../../Sonstiges/SimKernel/sim_main.h"

#ifndef INLINE
#define INLINE extern inline
#endif

#define REPORTS 25
#undef SHOW_DETAILED_PROGRESS

#define RESULT_DIMENSION 1

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
	unsigned int bins;
	unsigned int rawdatabins;
	// unsigned long mag der Kernel irgendwie nicht?
	long samples;
	unsigned int word_length;
	double std_noise;
	string inputfile_name;
	string outputfile_results_name;
	string outputfile_pars_name;
	gsl_rng* GSLrandom;
	double input_scaling;
	double cutoff;
	double tauF;
	bool OverrideRescalingQ;
	bool HighPassFilterQ;
	int lag;

	unsigned long *F_Ipast, *F_Jpast;
	unsigned long **F_Ipast_Jpast;

  rawdata **xdata;
#if RESULT_DIMENSION == 1
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
		sim.get("rawdatabins",rawdatabins);
		sim.get("bins",bins);
		// assert(bins==RESULT_DIMENSION);
		sim.get("samples",samples);
		sim.get("p",word_length,1,NoWarning);
		assert(word_length == 1);
		sim.get("noise",std_noise,0.0);
		sim.get("appliedscaling",input_scaling,1.0);
		sim.get("cutoff",cutoff,-1);
		sim.get("tauF",tauF,0);
		sim.get("OverrideRescalingQ",OverrideRescalingQ,false);
		sim.get("HighPassFilterQ",HighPassFilterQ,false);
		sim.get("lag",lag,0);
		assert(lag>=0);
		
		sim.get("inputfile",inputfile_name);
		sim.get("outputfile",outputfile_results_name);
		sim.get("outputparsfile",outputfile_pars_name);

		// initialize random number generator
		gsl_rng_env_setup();
		GSLrandom = gsl_rng_alloc(gsl_rng_default);
		gsl_rng_set(GSLrandom, 1234);
		
	};

	void execute(Sim& sim)
	{
	  cout <<"------ mutualinformation:main ------ olav, Mon 05 Apr 2010 ------"<<endl;
	  time_t start, end;
#ifdef SHOW_DETAILED_PROGRESS
	  time_t middle;
#endif
	  unsigned long totaltrials, completedtrials;

	  time(&start);
	  cout <<"start: "<<ctime(&start)<<flush;
	  cout <<"running on host: "<<flush;
	  system("hostname");
	  cout <<"current directory: "<<flush;
	  system("pwd");

	  cout <<"input file: "<<inputfile_name<<endl;
	  cout <<"output file: "<<outputfile_results_name<<endl;

	  cout <<"allocating memory..."<<flush;
	  xdata = new rawdata*[size];
#if RESULT_DIMENSION == 1
	  xresult = new double*[size];
#else
		xresult = new double**[size];
#endif
	  for(int i=0; i<size; i++)
	  {
	    xdata[i] = new rawdata[samples];
	    memset(xdata[i], 0, samples*sizeof(rawdata));
#if RESULT_DIMENSION == 1
	    xresult[i] = new double[size];
	    memset(xresult[i], 0, size*sizeof(double));
#else
			xresult[i] = new double*[size];
			for(int j=0; j<size; j++)
			{
				xresult[i][j] = new double[RESULT_DIMENSION];
				memset(xresult[i][j], 0, RESULT_DIMENSION*sizeof(double));
			}
#endif
	  }
		F_Ipast = new unsigned long[bins];
		F_Jpast = new unsigned long[bins];
		F_Ipast_Jpast = new unsigned long*[bins];
		for (rawdata x=0; x<bins; x++)
			F_Ipast_Jpast[x] = new unsigned long[bins];
		
	  cout <<" done."<<endl;

	  if (!OverrideRescalingQ)
			cout <<"loading data and adding noise (std "<<std_noise<<", cutoff "<<cutoff<<")..."<<flush;
		else cout <<"loading raw data (OverrideRescalingQ enabled)..."<<flush;
	  load_data();
	  cout <<" done."<<endl;

	  // main loop:
	  totaltrials = size*(size-1);
	  cout <<"set-up: "<<size<<" neurons, "<<samples<<" samples, "<<bins<<" bins, lag "<<lag<<endl;
		cout <<"assumed length of Markov chain: "<<word_length<<endl;
		// unsigned long long terms_sum = 0;
		// unsigned long long terms_zero = 0;
	
	  cout <<"running "<<flush;
	  for(int ii=0; ii<size; ii++)
	  {
	  	status(ii,REPORTS,size);
	    for(int jj=0; jj<size; jj++)
	    {
	      if (ii != jj)
	      {
	      	xresult[ii][jj] = MutualInformation(xdata[ii], xdata[jj]);
					// memset(tempresult, 0, bins*sizeof(double));
	      	// TransferEntropy(xdata[ii], xdata[jj],tempresult);
					// for (int k=0; k<bins; k++)
						// xresult[ii][jj][k] = tempresult[k];
	      }
	      else
	      {
#if RESULT_DIMENSION == 1
	      	xresult[ii][jj] = 0.0;
#endif
	      }
	    }
	  }
	  cout <<" done."<<endl;

	  time(&end);
	  cout <<"end: "<<ctime(&end)<<flush;
	  cout <<"runtime: "<<sec2string(difftime(end,start))<<endl;
	};
	
	void finalize(Sim& sim)
	{
#if RESULT_DIMENSION <= 1
	  write_result(xresult);
#else
		write_multidim_result(xresult);
#endif
		save_parameters();

		gsl_rng_free(GSLrandom);
				
		sim.io <<"End of Kernel (iteration="<<(sim.iteration())<<")"<<Endl;
	};
	

	double MutualInformation(rawdata *arrayI, rawdata *arrayJ)
	{
		/* see for reference:
		     Gourevitch und Eggermont. Evaluating Information Transfer Between Auditory
		     Cortical Neurons. Journal of Neurophysiology (2007) vol. 97 (3) pp. 2533 */
	  double result = 0.0;
		// memset(result, 0, bins*sizeof(double));

		// We are looking at the information flow of array1 ("I") -> array2 ("J")
			
		// clear memory
		memset(F_Ipast, 0, bins*sizeof(unsigned long));
		memset(F_Jpast, 0, bins*sizeof(unsigned long));
		for (rawdata x=0; x<bins; x++)
			memset(F_Ipast_Jpast[x], 0, bins*sizeof(unsigned long));
	
	  // extract probabilities (actually the number of occurrence)
		for (unsigned long t=lag; t<samples; t++)
	  {
			F_Ipast[arrayI[t-lag]]++;
			F_Jpast[arrayJ[t]]++;
			F_Ipast_Jpast[arrayI[t-lag]][arrayJ[t]]++;
		}
	
		for (rawdata k=0; k<bins; k++)
			for (rawdata l=0; l<bins; l++)
			{
				if (F_Ipast_Jpast[k][l] > 0)
					result += F_Ipast_Jpast[k][l]/double(samples-lag) * \
						log((samples-lag)*double(F_Ipast_Jpast[k][l])/(double(F_Ipast[k])*F_Jpast[l]));
			}
		
	  return result/log(2);
	};

	void load_data()
	{
		char* name = new char[inputfile_name.length()+1];
		strcpy(name,inputfile_name.c_str());
	  ifstream binaryfile(name, ios::binary);
		delete[] name;
		
		char* temparray = new char[samples];
		double* tempdoublearray = new double[samples];
		// memset(tempdoublearray, 0, samples*sizeof(double));
		int* tempintarray = new int[samples];

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

		// if(HighPassFilterQ)
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
			else
			// OverrideRescalingQ = false
			{
		    for(long k=0; k<samples; k++)
				{
					// transform to unsigned notation
					tempdoublearray[k] = double(temparray[k]);
					if (temparray[k]<0) tempdoublearray[k] += 256.;
					// transform back to original signal and apply noise (same as in Granger case)
					tempdoublearray[k] /= input_scaling;
					tempdoublearray[k] += gsl_ran_gaussian(GSLrandom,std_noise);
					
					// apply cutoff
					if ((cutoff>0)&&(tempdoublearray[k]>cutoff)) tempdoublearray[k] = cutoff;
				}
				
				if(HighPassFilterQ) {
					// of course, this is just a difference signal, so not really filtered
					memcpy(tempdoublearraycopy,tempdoublearray,samples*sizeof(double));
					tempdoublearray[0] = 0.0;
			    for(long k=1; k<samples; k++)
						tempdoublearray[k] = tempdoublearraycopy[k] - tempdoublearraycopy[k-1];
				}
					
				discretize(tempdoublearray,xdata[j]);
			}
	  }
	
		// cout <<endl;
		// for(int j=0; j<3000; j++)
		// 	cout <<"t = "<<double(j*tauF)<<" ms : xresponse[3] = "<<int(xdata[2][j])<<endl;
		// cout <<endl;
		// exit(1);

		delete[] temparray;
		delete[] tempdoublearray, tempdoublearraycopy;
	};

	void discretize(double* in, rawdata* out)
	{
		discretize(in,out,smallest(in,samples),largest(in,samples),bins);
	};
	void discretize(double* in, rawdata* out, unsigned int nr_bins)
	{
		discretize(in,out,smallest(in,samples),largest(in,samples),nr_bins);
	};
	void discretize(double* in, rawdata* out, double min, double max, unsigned int nr_bins)
	{
		// target binning is assumed to be 'bins'
		double xstepsize = (max-min)/(nr_bins-1);
		// cout <<"max = "<<max<<endl;
		// cout <<"min = "<<min<<endl;
		// cout <<"stepsize = "<<xstepsize<<endl;

		int xint;
		for (unsigned long t=0; t<samples; t++)
		{
			xint = round((in[t]-min)/xstepsize);
			// crop overshoot
			if (xint>=nr_bins) xint = bins-1;
			if (xint<0) xint = 0;

			out[t] = rawdata(xint);
			assert(out[t]<nr_bins);
		}
	};

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

		fileout1.precision(6);		
		fileout1 <<"{";
		fileout1 <<"executable->transferentropysim";
		fileout1 <<",iteration->"<<iteration;
		
		fileout1 <<",size->"<<size;
		fileout1 <<",rawdatabins->"<<rawdatabins;
		fileout1 <<",bins->"<<bins;
		fileout1 <<",cutoff->"<<cutoff;
		fileout1 <<",samples->"<<samples;
		fileout1 <<",p->"<<word_length;
		fileout1 <<",noise->"<<std_noise;
		fileout1 <<",tauF->"<<tauF;
		fileout1 <<",OverrideRescalingQ->"<<OverrideRescalingQ;
		fileout1 <<",HighPassFilterQ->"<<HighPassFilterQ;
		fileout1 <<",lag->"<<lag;
		
		fileout1 <<"}"<<endl;
		
		fileout1.close();
	};

	void write_result(double **array)
	{
		char* name = new char[outputfile_results_name.length()+1];
		strcpy(name,outputfile_results_name.c_str());
		ofstream fileout1(name);
		delete[] name;

		fileout1.precision(9);
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

	void write_multidim_result(double ***array)
	{
		char* name = new char[outputfile_results_name.length()+1];
		strcpy(name,outputfile_results_name.c_str());
		ofstream fileout1(name);
		delete[] name;

		fileout1.precision(9);
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
		    for(int k=0; k<RESULT_DIMENSION; k++)
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

	void generate_global(rawdata** raw, rawdata* global)
	{
		double avg;
		for (unsigned long t=0; t<samples; t++)
		{
			avg = 0.0;
			for (unsigned long j=0; j<size; j++)
				avg += double(raw[j][t]);
			global[t] = rawdata(round(avg/size));
			// test:
			// global[t] = 3;
		}
	};
};
