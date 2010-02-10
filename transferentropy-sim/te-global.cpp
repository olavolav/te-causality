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
#include "../../Simulationen/SimKernel/sim_main.h"

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
	unsigned int bins, globalbins;
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

	unsigned long** F_Ipast_Gpast;
	unsigned long*** F_Inow_Ipast_Gpast;
	unsigned long*** F_Ipast_Jpast_Gpast;
	unsigned long**** F_Inow_Ipast_Jpast_Gpast;
  rawdata **xdata;
	rawdata *xglobal;
  double **xresult;

  void initialize(Sim& sim)
	{
		iteration = sim.iteration();
		sim.io <<"Init: iteration "<<iteration<<", process "<< sim.process()<<Endl;
		time(&now);
		sim.io <<"ETA: "<<ETAstring(sim.iteration(),sim.n_iterations(),difftime(now,start))<<Endl;
		
		// read parameters from control file
		sim.get("size",size);
		sim.get("rawdatabins",rawdatabins);
		sim.get("bins",bins);
		sim.get("globalbins",globalbins);
		sim.get("samples",samples);
		sim.get("p",word_length);
		assert(word_length == 1);
		sim.get("noise",std_noise);
		sim.get("appliedscaling",input_scaling);
		
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
	  cout <<"------ transferentropy-sim:global ------ olav, Wed 10 Jun 2009 ------"<<endl;
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
	
		F_Ipast_Gpast = new unsigned long*[bins];
		F_Inow_Ipast_Gpast = new unsigned long**[bins];
		F_Ipast_Jpast_Gpast = new unsigned long**[bins];
		F_Inow_Ipast_Jpast_Gpast = new unsigned long***[bins];
		for (char x=0; x<bins; x++)
		{
			F_Ipast_Gpast[x] = new unsigned long[globalbins];
			F_Inow_Ipast_Gpast[x] = new unsigned long*[bins];
			F_Ipast_Jpast_Gpast[x] = new unsigned long*[bins];
			F_Inow_Ipast_Jpast_Gpast[x] = new unsigned long**[bins];
			for (char x2=0; x2<bins; x2++)
			{
				F_Inow_Ipast_Gpast[x][x2] = new unsigned long[globalbins];
				F_Ipast_Jpast_Gpast[x][x2] = new unsigned long[globalbins];
				F_Inow_Ipast_Jpast_Gpast[x][x2] = new unsigned long*[bins];
				for (char x3=0; x3<bins; x3++)
					F_Inow_Ipast_Jpast_Gpast[x][x2][x3] = new unsigned long[globalbins];
			}
		}	
	  cout <<" done."<<endl;

	  cout <<"loading data, adding noise (std "<<std_noise<<") and generating global signal... "<<flush;
		xglobal = new rawdata[samples];
	  load_data();
	  // generate_global();
	  cout <<" done."<<endl;
	
	  // main loop:
	  totaltrials = size*(size-1);
	  cout <<"set-up: "<<size<<" neurons, "<<samples<<" samples, "<<bins<<" bins, "<<globalbins<<" globalbins"<<endl;
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
		for (char x=0; x<bins; x++)
		{
			memset(F_Ipast_Gpast[x], 0, globalbins*sizeof(unsigned long));
			for (char x2=0; x2<bins; x2++)
			{
				memset(F_Inow_Ipast_Gpast[x][x2], 0, globalbins*sizeof(unsigned long));
				memset(F_Ipast_Jpast_Gpast[x][x2], 0, globalbins*sizeof(unsigned long));
				for (char x3=0; x3<bins; x3++)
					memset(F_Inow_Ipast_Jpast_Gpast[x][x2][x3], 0, globalbins*sizeof(unsigned long));
			}
		}
	
	  // extract probabilities (actually number of occurrence)
		for (unsigned long t=word_length; t<samples; t++)
	  {
			F_Ipast_Gpast[arrayI[t-1]][xglobal[t-1]]++;
			F_Inow_Ipast_Gpast[arrayI[t]][arrayI[t-1]][xglobal[t-1]]++;
			F_Ipast_Jpast_Gpast[arrayI[t-1]][arrayJ[t-1]][xglobal[t-1]]++;
			F_Inow_Ipast_Jpast_Gpast[arrayI[t]][arrayI[t-1]][arrayJ[t-1]][xglobal[t-1]]++;
	#ifdef SHOW_DETAILED_PROGRESS
			status(t, REPORTS, samples-word_length);
	#endif
		}

		// index convention:
		// k - Ipast
		// l - Jpast
		// m - Inow
		// g - Gpast

		// calculate transfer entropy
		// for (char k=0; k<bins; k++)
		// 	for (char g=0; g<bins; g++)
		// 		if (F_Ipast_Gpast[k][g]!=0) for (char l=0; l<bins; l++)
		// 			if (F_Ipast_Jpast_Gpast[k][l][g]!=0) for (char m=0; m<bins; m++)
		// 					if ((F_Inow_Ipast_Gpast[m][k][g]!=0) && (F_Inow_Ipast_Jpast_Gpast[m][k][l][g] != 0))
		// 						result += double(F_Inow_Ipast_Jpast_Gpast[m][k][l][g])/double(samples-word_length) * log(double(F_Inow_Ipast_Jpast_Gpast[m][k][l][g])*double(F_Ipast_Gpast[k][g])/(double(F_Ipast_Jpast_Gpast[k][l][g])*double(F_Inow_Ipast_Gpast[m][k][g])));
		double Hxx = 0.0;
		for (char k=0; k<bins; k++)
			for (char m=0; m<bins; m++)
				for (char g=0; g<globalbins; g++)
					if (F_Inow_Ipast_Gpast[m][k][g] > 0)
						Hxx -= double(F_Inow_Ipast_Gpast[m][k][g])/(samples-word_length) * log(double(F_Inow_Ipast_Gpast[m][k][g])/double(F_Ipast_Gpast[k][g]));
		Hxx /= log(2);

		double Hxxy = 0.0;
		for (char k=0; k<bins; k++)
			for (char l=0; l<bins; l++)
				for (char m=0; m<bins; m++)
					for (char g=0; g<globalbins; g++)
						if (F_Inow_Ipast_Jpast_Gpast[m][k][l][g] > 0)
							Hxxy -= double(F_Inow_Ipast_Jpast_Gpast[m][k][l][g])/(samples-word_length) * log(double(F_Inow_Ipast_Jpast_Gpast[m][k][l][g])/double(F_Ipast_Jpast_Gpast[k][l][g]));
		Hxxy /= log(2);

	  return (Hxx - Hxxy);
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

	  for(int j=0; j<size; j++)
	  {
	    binaryfile.read(temparray, samples);
	    for(long k=0; k<samples; k++)
			{
				// transform to unsigned notation
				tempdoublearray[k] = double(temparray[k]);
				if (temparray[k]<0) tempdoublearray[k] += 256.;
				// transform back to original signal and apply noise (same as in Granger case)
				tempdoublearray[k] /= input_scaling;
				tempdoublearray[k] += gsl_ran_gaussian(GSLrandom,std_noise);
				
				// rescaling to new bins (because rawdatabins/input_scaling is the old maximum)
				// xtemp = xtemp/(rawdatabins/input_scaling)*bins;
				// xtemp = xtemp/0.45*bins;
				// if (xtemp < 0) xtemp = 0.0;
				// if (xtemp > bins-1) xtemp = bins-1;
				
				// add to what later becomes the global signal
				xglobaltemp[k] += tempdoublearray[k];
				
				// convert to target data format
				// xdata[j][k] = discretize(xtemp, /* 3/4*rawdatabins/input_scaling,*/ 2*std_noise);
				// xdata[j][k] = rawdata(round(xtemp));
			}
			discretize(tempdoublearray,xdata[j]);
	  }
	
		// cout <<endl;
		// for(int j=0; j<400; j++)
		// 	cout <<int(xdata[2][j])<<",";
		// cout <<endl;
		// exit(1);

		// generate global signal (rescaled to globalbins binning)
		for (unsigned long t=0; t<samples; t++)
		{
			xglobaltemp[t] /= size;
			// xglobal[t] = rawdata(round(xglobaltemp[t]/size*globalbins/bins));
			// xglobal[t] = 0; // for testing
			// assert(xglobal[t]<globalbins);
			// if (xglobal[t] >= globalbins) xglobal[t] = globalbins-1;
		}
		discretize(xglobaltemp,xglobal,globalbins);
	
		delete[] temparray;
		delete[] xglobaltemp;
		delete[] tempdoublearray;
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
		fileout1 <<"executable->teglobalsim";
		fileout1 <<",iteration->"<<iteration;
		
		fileout1 <<",size->"<<size;
		fileout1 <<",rawdatabins->"<<rawdatabins;
		fileout1 <<",bins->"<<bins;
		fileout1 <<",globalbins->"<<globalbins;
		fileout1 <<",samples->"<<samples;
		fileout1 <<",p->"<<word_length;
		fileout1 <<",noise->"<<std_noise;
		
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

	/* void generate_global()
	{
		double avg;
		for (unsigned long t=0; t<samples; t++)
		{
			avg = 0.0;
			for (unsigned long j=0; j<size; j++)
				avg += double(xdata[j][t]);
			xglobal[t] = rawdata(floor(avg/size*globalbins/bins));
			assert(xglobal[t]<globalbins);
			// test:
			// global[t] = 3;
		}
	}; */
};
