// calculate the Granger causality between a numer of time series
// created by olav, Mi  2 Dez 2009 21:52:52 CET

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <ctime>
#include <cstring>
#include <sstream>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "../../Simulationen/olav.h"
#include "../../Simulationen/SimKernel/sim_main.h"

#define REPORTS 25
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
	// unsigned long mag der Kernel irgendwie nicht?
	long samples;
	unsigned int regression_order;
	double std_noise;
	string inputfile_name;
	string outputfile_results_name;
	string outputfile_pars_name;
	gsl_rng* GSLrandom;
	double input_scaling;
	double cutoff;
	double tauF;
	
	double** xdata;
	double* xglobal;
	int detrend_mode;
#if RESULT_DIMENSION > 1
	double ***xresult;
#else
	double **xresult;
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
		sim.get("bins",bins);
		sim.get("samples",samples);
		sim.get("p",regression_order);
		sim.get("noise",std_noise);
		sim.get("appliedscaling",input_scaling);
		sim.get("detrend",detrend_mode);
		assert((detrend_mode==0)||(detrend_mode==1));

		sim.get("cutoff",cutoff);
		sim.get("tauF",tauF,0);
		
		sim.get("inputfile",inputfile_name);
		sim.get("outputfile",outputfile_results_name);
		sim.get("outputparsfile",outputfile_pars_name);
		
		// initialize random number generator
		gsl_rng_env_setup();
		GSLrandom = gsl_rng_alloc(gsl_rng_default);
		gsl_rng_set(GSLrandom, 1234);
		
		// cout <<"bins = "<<bins<<", samples= "<<samples<<endl;
	};
	
	void execute(Sim& sim)
	{
	  cout <<"------ granger-sim:main ------ olav, Fri 29 Jan 2010 ------"<<endl;
	  time_t start, end;
	  time_t middle;

	  time(&start);
	  cout <<"start: "<<ctime(&start)<<flush;
	  cout <<"running on host: "<<flush;
	  system("hostname");
	  cout <<"current directory: "<<flush;
	  system("pwd");

	  cout <<"input file: "<<inputfile_name<<endl;
	  cout <<"output file: "<<outputfile_results_name<<endl;
	  cout <<"pars output file: "<<outputfile_pars_name<<endl;

	  cout <<"allocating memory..."<<flush;
	  xdata = new double*[size];
	  for(int i=0; i<size; i++)
	    xdata[i] = new double[samples];
		if (detrend_mode == 1)
			xglobal = new double[samples];

#if RESULT_DIMENSION > 1
	  xresult = new double**[size];
	  for(int i=0; i<size; i++)
		{
	    xresult[i] = new double*[size];
	    for(int j=0; j<size; j++)
				xresult[i][j] = new double[regression_order+1];
	  }
#else
	  xresult = new double*[size];
	  for(int i=0; i<size; i++)
	    xresult[i] = new double[size];
#endif

		// allocate GSL workspaces
		gsl_multifit_linear_workspace * GSLworkspaceBoth = \
			gsl_multifit_linear_alloc(samples-regression_order,2*regression_order+1);
		double residue;
		// the coulmns in input are: (nr. of columns)
		// target (regression_order), source (regression_order), constant term (1)
		gsl_matrix* inputBoth = gsl_matrix_alloc(samples-regression_order,2*regression_order+1);
		// because the last column of this matrix never changes, we can set it here:
		gsl_matrix_set_all(inputBoth,1.0);
		gsl_matrix* covBoth = gsl_matrix_alloc(2*regression_order+1,2*regression_order+1);
		// output is the only GSL consturuct that we use for both cases ("both" and "single")
		gsl_vector* output = gsl_vector_alloc(samples-regression_order);
		gsl_vector* coeffBoth = gsl_vector_alloc(2*regression_order+1);

		// gsl_multifit_linear_workspace * GSLworkspaceSingle =	gsl_multifit_linear_alloc(samples-regression_order,1*regression_order+1);
		// gsl_matrix* inputSingle = gsl_matrix_alloc(samples-regression_order,1*regression_order+1);
		// gsl_matrix_set_all(inputSingle,1.0);
		// gsl_matrix* covSingle = gsl_matrix_alloc(1*regression_order+1,1*regression_order+1);
		// gsl_vector* coeffSingle = gsl_vector_alloc(1*regression_order+1);

	  cout <<" done."<<endl;

	  cout <<"loading data and adding noise (std "<<std_noise<<", cutoff "<<cutoff<<")..."<<flush;
	  load_data();
	  cout <<" done."<<endl;
	
		if (detrend_mode==1)
		{
			cout <<"generating global signal..."<<flush;
		  generate_global();
		  cout <<" done."<<endl;
		
			cout <<"detrending time series..."<<flush;
		  detrend_signals();
		  cout <<" done."<<endl;		  
		}

		cout <<"set-up: "<<size<<" neurons";
		cout <<", "<<size*(size-1)<<" trials";
		cout <<", regression order "<<regression_order;
		cout <<", "<<samples<<" samples"<<endl;

		// double test_residue;
		// double covTraceBoth, covTraceSingle;
	  // main loop:
	  for(int ii=0; ii<size; ii++)
	  {
			cout <<"node #"<<ii+1<<": "<<flush;

			// set up source data of connection
			for(long tt=0; tt<samples-regression_order; tt++)
				for(long tt2=0; tt2<regression_order; tt2++)
				{
					gsl_matrix_set(inputBoth,tt,regression_order+tt2,xdata[ii][tt-tt2]);
					// gsl_matrix_set(inputSingle,tt,tt2,xdata[ii][tt-tt2]);
				}

	    for(int jj=0; jj<size; jj++)
	    {
		  	status(jj,REPORTS,size);
	      if (ii != jj)
	      {
					// add target data of connection and perform auto-regression
					for(long tt=0; tt<samples-regression_order; tt++)
					{
						gsl_vector_set(output,tt,xdata[jj][tt+1]); // is this correct for orders != 1?

						for(long tt2=0; tt2<regression_order; tt2++)
							gsl_matrix_set(inputBoth,tt,tt2,xdata[jj][tt-tt2]);
					}

					gsl_multifit_linear(inputBoth,output,coeffBoth,covBoth,&residue,GSLworkspaceBoth);
					// test_residue = residue;
					// gsl_multifit_linear(inputSingle,output,coeffSingle,covSingle,&residue,GSLworkspaceSingle);

					// test:
					// xresult[ii][jj] = residue/test_residue;

					// save those coefficients that depend on the source
					// covTraceBoth = covTraceSingle = 0.0;
					// for (int kk=0; kk<regression_order; kk++)
					// {
						// xresult[ii][jj][kk] = gsl_vector_get(coeffBoth,regression_order+kk);
						// covTraceBoth += gsl_matrix_get(covBoth,regression_order+kk,regression_order+kk);
						// covTraceSingle += gsl_matrix_get(covSingle,regression_order+kk,regression_order+kk);
					// }
					// xresult[ii][jj][regression_order] = residue;
					xresult[ii][jj] = gsl_vector_get(coeffBoth,regression_order);

					// save reduction in variance as xresult
					// xresult[ii][jj] = log(sqrt(gsl_matrix_get(covSingle,0,0)+gsl_matrix_get(covSingle,1,1))/sqrt(gsl_matrix_get(covBoth,0,0)+gsl_matrix_get(covBoth,1,1)));
	      }
	      // else for (int kk=0; kk<regression_order; kk++)
					// xresult[ii][jj][kk] = 0.0;
				else xresult[ii][jj] = 0.0;
				// cout <<"debug: result="<<xresult[ii][jj]<<endl;
	    }
			time(&middle);
			cout <<" (elapsed "<<sec2string(difftime(middle,start))<<",";
			cout <<" ETA "<<ETAstring(ii+1,size,difftime(middle,start))<<")"<<endl;
	  }

	  time(&end);
	  cout <<"end: "<<ctime(&end)<<flush;
	  cout <<"runtime: "<<sec2string(difftime(end,start))<<endl;

		save_parameters();
#if RESULT_DIMENSION > 1
	  write_multidim_result(xresult);
#else
	  write_result(xresult);
#endif

	  // free memory
		gsl_multifit_linear_free(GSLworkspaceBoth);
		gsl_matrix_free(inputBoth);
		gsl_matrix_free(covBoth);
		gsl_vector_free(output);
		gsl_vector_free(coeffBoth);
		// gsl_multifit_linear_free(GSLworkspaceSingle);
		// gsl_matrix_free(inputSingle);
		// gsl_matrix_free(covSingle);
		// gsl_vector_free(coeffSingle);
		for(int i=0; i<size; i++) delete[] xdata[i];
	  delete[] xdata;
	  for(int i=0; i<size; i++) delete[] xresult[i];
	  delete[] xresult;		
	};

	void finalize(Sim& sim)
	{
		gsl_rng_free(GSLrandom);
		
		sim.io <<"End of Kernel (iteration="<<(sim.iteration())<<")"<<Endl;
	};
	
	void load_data()
	{
		char* name = new char[inputfile_name.length()+1];
		strcpy(name,inputfile_name.c_str());
	  ifstream binaryfile(name, ios::binary);
		delete[] name;
		
		char* temparray = new char[samples];
		double xtemp;

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
				xtemp = double(temparray[k]);
				if (temparray[k]<0) xtemp += 256.;
				// apply noise (same as in Granger case)
				xdata[j][k] = xtemp/input_scaling + gsl_ran_gaussian(GSLrandom,std_noise);
				// apply cutoff
				if ((cutoff>0)&&(xdata[j][k]>cutoff)) xdata[j][k] = cutoff;
			}
	  }
	
		// for(int t=700;t<900;t++)
		// 	cout <<"frame "<<t<<": "<<xdata[3][t]<<endl;
		// exit(1);
	
		delete[] temparray;
	};
	
	void save_parameters()
	{
		char* name = new char[outputfile_pars_name.length()+1];
		strcpy(name,outputfile_pars_name.c_str());
		ofstream fileout1(name);
		delete[] name;

		fileout1.precision(6);		
		fileout1 <<"{";
		fileout1 <<"iteration->"<<iteration;
		
		fileout1 <<",detrend->"<<detrend_mode;
		fileout1 <<",size->"<<size;
		fileout1 <<",bins->"<<bins;
		fileout1 <<",samples->"<<samples;
		fileout1 <<",inputscaling->"<<input_scaling;
		fileout1 <<",cutoff->"<<cutoff;
		fileout1 <<",p->"<<regression_order;
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
	    	fileout1 <<array[j][i];
	    }
	    fileout1 <<"}"<<endl;
	  }
	  fileout1 <<"}"<<endl;

	  cout <<"Granger causality matrix saved."<<endl;
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

	  cout <<"Granger causality matrix saved."<<endl;
	};
	
	void generate_global()
	{
		double avg;
		for (unsigned long t=0; t<samples; t++)
		{
			avg = 0.0;
			for (unsigned long j=0; j<size; j++)
				avg += double(xdata[j][t]);
			xglobal[t] = avg/size;
		}
	};
	
	void detrend_signals()
	{
		for (unsigned long t=0; t<samples; t++)
			for (unsigned long j=0; j<size; j++)
				xdata[j][t] -= xglobal[t];
	};
	
};
