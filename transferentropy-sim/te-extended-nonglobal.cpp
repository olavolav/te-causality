// calculate the transfer entropy between a numer of time series
// created by olav, Mi 10 Jun 2009 19:42:11 IDT

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <ctime>
#include <cstring>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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

	unsigned long* F_Ipast;
	unsigned long** F_Inow_Ipast;
	unsigned long*** F_Ipast_Jpast_Jpast2;
	unsigned long**** F_Inow_Ipast_Jpast_Jpast2;
  rawdata **xdata;
  double **xresult;

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
		sim.get("samples",samples);
		sim.get("p",word_length);
		assert(word_length == 1);
		sim.get("noise",std_noise);
		sim.get("appliedscaling",input_scaling);
		sim.get("cutoff",cutoff);
		sim.get("tauF",tauF);
		
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
	  cout <<"------ transferentropy-sim:extended ------ olav, Mon 15 Feb 2010 ------"<<endl;
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
	  xresult = new double*[size];
	  for(int i=0; i<size; i++)
	  {
	    xdata[i] = new rawdata[samples];
	    memset(xdata[i], 0, samples*sizeof(rawdata));
	    xresult[i] = new double[size];
	    memset(xresult[i], 0, size*sizeof(double));
	  }
		F_Ipast = new unsigned long[bins];
		F_Inow_Ipast = new unsigned long*[bins];
		F_Ipast_Jpast_Jpast2 = new unsigned long**[bins];
		F_Inow_Ipast_Jpast_Jpast2 = new unsigned long***[bins];
		for (rawdata x=0; x<bins; x++)
		{
			F_Inow_Ipast[x] = new unsigned long[bins];

			F_Ipast_Jpast_Jpast2[x] = new unsigned long*[bins];
			F_Inow_Ipast_Jpast_Jpast2[x] = new unsigned long**[bins];
			for (rawdata x2=0; x2<bins; x2++)
			{
				F_Ipast_Jpast_Jpast2[x][x2] = new unsigned long[bins];
				F_Inow_Ipast_Jpast_Jpast2[x][x2] = new unsigned long*[bins];
				for (rawdata x3=0; x3<bins; x3++)
				{
					F_Inow_Ipast_Jpast_Jpast2[x][x2][x3] = new unsigned long[bins];
				}
			}
		}
	  cout <<" done."<<endl;

	  cout <<"loading data and adding noise (std "<<std_noise<<", cutoff "<<cutoff<<")..."<<flush;
	  load_data();
	  cout <<" done."<<endl;

	  // main loop:
	  totaltrials = size*(size-1);
	  cout <<"set-up: "<<size<<" neurons, "<<samples<<" samples, "<<bins<<" bins"<<endl;
		cout <<"assumed length of Markov chain: "<<word_length<<endl;
	  completedtrials = 0;
		// unsigned long long terms_sum = 0;
		// unsigned long long terms_zero = 0;
	
#ifndef SHOW_DETAILED_PROGRESS
	  	cout <<"running "<<flush;
#endif

	  for(int ii=0; ii<size; ii++)
	  {
#ifndef SHOW_DETAILED_PROGRESS
	  	status(ii,REPORTS,size);
#endif
	    for(int jj=0; jj<size; jj++)
	    {
	      if (ii != jj)
	      {
#ifdef SHOW_DETAILED_PROGRESS
	      	cout <<"#"<<ii+1<<" -> #"<<jj+1<<": "<<flush;
#endif
	      	xresult[ii][jj] = TransferEntropy(xdata[ii], xdata[jj]);
					completedtrials++;
#ifdef SHOW_DETAILED_PROGRESS
					time(&middle);
					cout <<" (elapsed "<<sec2string(difftime(middle,start))<<",";
					// cout <<" ETA "<<sec2string((totaltrials-(completedtrials+1))*difftime(middle,start)/(completedtrials+1))<<")"<<endl;
					cout <<" ETA "<<ETAstring(completedtrials,totaltrials,difftime(middle,start))<<")"<<endl;
					// cout <<" (elapsed: "<<sec2string(difftime(middle,start))<<")"<<endl;
#endif
	      }
	      else
	      {
#ifdef SHOW_DETAILED_PROGRESS
					cout <<"#"<<ii+1<<" -> #"<<jj+1<<": skipped."<<endl;
#endif
	      	xresult[ii][jj] = 0.0;
	      }
	    }
	  }
#ifndef SHOW_DETAILED_PROGRESS
	  cout <<" done."<<endl;
#endif

	  time(&end);
	  cout <<"end: "<<ctime(&end)<<flush;
	  cout <<"runtime: "<<sec2string(difftime(end,start))<<endl;

		// cout <<"TE terms: "<<terms_sum<<", of those zero: "<<terms_zero<<" ("<<int(double(terms_zero)*100/terms_sum)<<"%)"<<endl;
	};
	
	void finalize(Sim& sim)
	{
	  write_result(xresult);
		save_parameters();

		gsl_rng_free(GSLrandom);
				
		sim.io <<"End of Kernel (iteration="<<(sim.iteration())<<")"<<Endl;
	};
	

	double TransferEntropy(rawdata *arrayI, rawdata *arrayJ)
	{
		/* see for reference:
		     Gourevitch und Eggermont. Evaluating Information Transfer Between Auditory
		     Cortical Neurons. Journal of Neurophysiology (2007) vol. 97 (3) pp. 2533 */
	  double result = 0.0;

		// We are looking at the information flow of array1 ("I") -> array2 ("J")
	
		// allocate memory
		memset(F_Ipast, 0, bins*sizeof(unsigned long));
		for (rawdata x=0; x<bins; x++)
		{
			memset(F_Inow_Ipast[x], 0, bins*sizeof(unsigned long));
			for (rawdata x2=0; x2<bins; x2++)
			{
				memset(F_Ipast_Jpast_Jpast2[x][x2], 0, bins*sizeof(unsigned long));
				for (rawdata x3=0; x3<bins; x3++)
					memset(F_Inow_Ipast_Jpast_Jpast2[x][x2][x3], 0, bins*sizeof(unsigned long));
			}
		}
	
	  // extract probabilities (actually number of occurrence)
		for (unsigned long t=2; t<samples; t++)
	  {
			F_Ipast[arrayI[t-1]]++;
			F_Inow_Ipast[arrayI[t]][arrayI[t-1]]++;
			F_Ipast_Jpast_Jpast2[arrayI[t-1]][arrayJ[t-1]][arrayJ[t-2]]++;
			F_Inow_Ipast_Jpast_Jpast2[arrayI[t]][arrayI[t-1]][arrayJ[t-1]][arrayJ[t-2]]++;
#ifdef SHOW_DETAILED_PROGRESS
			status(t, REPORTS, samples-word_length);
#endif
		}
	
		for (rawdata k=0; k<bins; k++)
			for (rawdata l=0; l<bins; l++)
				for (rawdata m=0; m<bins; m++)
					for (rawdata l2=0; l2<bins; l2++)
						if (F_Inow_Ipast_Jpast_Jpast2[m][k][l][l2] != 0)
							result += F_Inow_Ipast_Jpast_Jpast2[m][k][l][l2]/double(samples-2) * \
								log(double(F_Inow_Ipast_Jpast_Jpast2[m][k][l][l2]*F_Ipast[k])/(F_Ipast_Jpast_Jpast2[k][l][l2]*F_Inow_Ipast[m][k]));
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
					
				if(!AdaptiveBinningQ) discretize(tempdoublearray,xdata[j],bins);
				else discretize2accordingtoStd(tempdoublearray,xdata[j]);
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
		fileout1 <<"executable->transferentropysim-extended";
		fileout1 <<",iteration->"<<iteration;
		
		fileout1 <<",size->"<<size;
		fileout1 <<",rawdatabins->"<<rawdatabins;
		fileout1 <<",bins->"<<bins;
		fileout1 <<",cutoff->"<<cutoff;
		fileout1 <<",samples->"<<samples;
		fileout1 <<",p->"<<word_length;
		fileout1 <<",noise->"<<std_noise;
		fileout1 <<",tauF->"<<tauF;
		
		fileout1 <<"}"<<endl;
		
		fileout1.close();
	};

	void write_result(double **array)
	{
		char* name = new char[outputfile_results_name.length()+1];
		strcpy(name,outputfile_results_name.c_str());
		ofstream fileout1(name);
		delete[] name;

		fileout1.precision(6);
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
