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
	gsl_permutation* GSLpermutation;
	double input_scaling;
	double cutoff;
	double tauF;
	bool OverrideRescalingQ;
	bool HighPassFilterQ;
	bool InstantFeedbackTermQ;
	bool GourevitchNormalizationQ;

	unsigned long* F_Ipast;
	unsigned long** F_Inow_Ipast;
	unsigned long** F_Ipast_Jpast;
	unsigned long*** F_Inow_Ipast_Jpast;
  rawdata **xdata;
	rawdata *xdatashuffel;
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
		sim.get("InstantFeedbackTermQ",InstantFeedbackTermQ,false);
		sim.get("GourevitchNormalizationQ",GourevitchNormalizationQ,false);
		
		sim.get("inputfile",inputfile_name);
		sim.get("outputfile",outputfile_results_name);
		sim.get("outputparsfile",outputfile_pars_name);

		// initialize random number generator
		gsl_rng_env_setup();
		GSLrandom = gsl_rng_alloc(gsl_rng_default);
		gsl_rng_set(GSLrandom, 1234);
		
		// initialize permutation generator
		if (GourevitchNormalizationQ)
			GSLpermutation = gsl_permutation_alloc(samples);
		
	};

	void execute(Sim& sim)
	{
	  cout <<"------ transferentropy-sim:main ------ olav, Wed 10 Jun 2009 ------"<<endl;
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
		if(GourevitchNormalizationQ)
			xdatashuffel = new rawdata[samples];
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
		F_Inow_Ipast = new unsigned long*[bins];
		F_Ipast_Jpast = new unsigned long*[bins];
		F_Inow_Ipast_Jpast = new unsigned long**[bins];
		for (rawdata x=0; x<bins; x++)
		{
			F_Inow_Ipast[x] = new unsigned long[bins];
			F_Ipast_Jpast[x] = new unsigned long[bins];

			F_Inow_Ipast_Jpast[x] = new unsigned long*[bins];
			for (rawdata x2=0; x2<bins; x2++)
			{
				F_Inow_Ipast_Jpast[x][x2] = new unsigned long[bins];
			}
		}
	  cout <<" done."<<endl;

	  if (!OverrideRescalingQ)
			cout <<"loading data and adding noise (std "<<std_noise<<", cutoff "<<cutoff<<")..."<<flush;
		else cout <<"loading raw data (OverrideRescalingQ enabled)..."<<flush;
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

		// double* tempresult = new double[bins];
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
	      	if(GourevitchNormalizationQ)
						xresult[ii][jj] = NormalizedTransferEntropy(xdata[ii], xdata[jj]);
					else xresult[ii][jj] = TransferEntropy(xdata[ii], xdata[jj]);
					// memset(tempresult, 0, bins*sizeof(double));
	      	// TransferEntropy(xdata[ii], xdata[jj],tempresult);
					// for (int k=0; k<bins; k++)
						// xresult[ii][jj][k] = tempresult[k];
	
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
#if RESULT_DIMENSION == 1
	      	xresult[ii][jj] = 0.0;
#endif
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
#if RESULT_DIMENSION <= 1
	  write_result(xresult);
#else
		write_multidim_result(xresult);
#endif
		save_parameters();

		// free allocated memory
		gsl_rng_free(GSLrandom);
		if (GourevitchNormalizationQ)
			gsl_permutation_free(GSLpermutation);
		
#if RESULT_DIMENSION > 1
		for (int x=0; x<RESULT_DIMENSION; x++)
			delete[] xresult[x];
#endif
		delete[] xresult;
		
		delete[] F_Ipast;
		for (rawdata x=0; x<bins; x++)
		{
			delete[] F_Inow_Ipast[x];
			delete[] F_Ipast_Jpast[x];
			for (rawdata x2=0; x2<bins; x2++)
				delete[] F_Inow_Ipast_Jpast[x][x2];
			delete[] F_Inow_Ipast_Jpast[x];
			delete[] xdata[x];
		}
		delete[] F_Inow_Ipast;
		delete[] F_Ipast_Jpast;
		delete[] F_Inow_Ipast_Jpast;
		delete[] xdata;
		
		sim.io <<"End of Kernel (iteration="<<(sim.iteration())<<")"<<Endl;
	};
	

	double TransferEntropy(rawdata *arrayI, rawdata *arrayJ)
	{
		/* see for reference:
		     Gourevitch und Eggermont. Evaluating Information Transfer Between Auditory
		     Cortical Neurons. Journal of Neurophysiology (2007) vol. 97 (3) pp. 2533 */
	  double result = 0.0;
		// memset(result, 0, bins*sizeof(double));

		// We are looking at the information flow of array1 ("I") -> array2 ("J")
			
		// clear memory
		memset(F_Ipast, 0, bins*sizeof(unsigned long));
		for (rawdata x=0; x<bins; x++)
		{
			memset(F_Inow_Ipast[x], 0, bins*sizeof(unsigned long));
			memset(F_Ipast_Jpast[x], 0, bins*sizeof(unsigned long));
			for (rawdata x2=0; x2<bins; x2++)
				memset(F_Inow_Ipast_Jpast[x][x2], 0, bins*sizeof(unsigned long));
		}
	
	  // extract probabilities (actually the number of occurrence)
		unsigned long const JShift = 0 + 1*InstantFeedbackTermQ;
		for (unsigned long t=word_length; t<samples; t++)
	  {
			F_Ipast[arrayI[t-1]]++;
			F_Inow_Ipast[arrayI[t]][arrayI[t-1]]++;
			F_Ipast_Jpast[arrayI[t-1]][arrayJ[t-1+JShift]]++;
			F_Inow_Ipast_Jpast[arrayI[t]][arrayI[t-1]][arrayJ[t-1+JShift]]++;
#ifdef SHOW_DETAILED_PROGRESS
			status(t, REPORTS, samples-word_length);
#endif
		}
	
		for (rawdata k=0; k<bins; k++)
			for (rawdata l=0; l<bins; l++)
			{
				if (F_Ipast_Jpast[k][l] > 0)
					for (rawdata m=0; m<bins; m++)
					// test: for (rawdata m=k+1; m<bins; m++)
						// if (F_Ipast[m]*F_Inow_Ipast[m][l]*F_Inow_Ipast_Jpast[m][k][l] != 0)
						if (F_Inow_Ipast_Jpast[m][k][l] != 0)
						{
							result += F_Inow_Ipast_Jpast[m][k][l]/double(samples-word_length) * \
								log(double(F_Inow_Ipast_Jpast[m][k][l]*F_Ipast[k])/(F_Ipast_Jpast[k][l]*F_Inow_Ipast[m][k]));
							// result[m] += F_Inow_Ipast_Jpast[m][k][l]/double(samples-word_length) * \
								log(double(F_Inow_Ipast_Jpast[m][k][l]*F_Ipast[k])/(F_Ipast_Jpast[k][l]*F_Inow_Ipast[m][k]));
						}
			}
		
		// for (rawdata k=0; k<bins; k++) result[k] /= log(2);		
	  return result/log(2);
	};
	
	double SingleSeriesTransitionEntropy(rawdata *arrayI)
	{
	  double result = 0.0;

		// clear memory
		memset(F_Ipast, 0, bins*sizeof(unsigned long));
		for (rawdata x=0; x<bins; x++)
			memset(F_Inow_Ipast[x], 0, bins*sizeof(unsigned long));

	  // extract probabilities (actually the number of occurrence)
		for (unsigned long t=word_length; t<samples; t++)
	  {
			F_Ipast[arrayI[t-1]]++;
			F_Inow_Ipast[arrayI[t]][arrayI[t-1]]++;
		}

		for (rawdata k=0; k<bins; k++)
			for (rawdata m=0; m<bins; m++)
				if (F_Inow_Ipast[m][k] != 0)
				{
					result += F_Inow_Ipast[m][k]/double(samples-word_length) * \
						log(double(F_Inow_Ipast[m][k])/F_Ipast[k]);
				}

	  return result/log(2);
	};
	

	double NormalizedTransferEntropy(rawdata *arrayI, rawdata *arrayJ)
	{
		/* see for reference:
		     Gourevitch und Eggermont. Evaluating Information Transfer Between Auditory
		     Cortical Neurons. Journal of Neurophysiology (2007) vol. 97 (3) pp. 2533 */
		double unnormalizedTE = TransferEntropy(arrayI,arrayJ);
		
		memcpy(xdatashuffel,arrayI,samples*sizeof(rawdata));
		gsl_ran_shuffle(GSLrandom,xdatashuffel,samples,sizeof(rawdata));
		double shuffeledTE = TransferEntropy(xdatashuffel,arrayJ);
		
		// inefficient
		double ownpastH = SingleSeriesTransitionEntropy(arrayJ);
		// double ownpastH = 1.0;
		
		return (unnormalizedTE-shuffeledTE)/ownpastH;
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
		fileout1 <<",InstantFeedbackTermQ->"<<InstantFeedbackTermQ;
		fileout1 <<",GourevitchNormalizationQ->"<<GourevitchNormalizationQ;
		
		fileout1 <<",inputfile->\""<<inputfile_name<<"\"";
		fileout1 <<",outputfile->\""<<outputfile_results_name<<"\"";
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
