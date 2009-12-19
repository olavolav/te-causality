// calculate the Granger causality between a numer of time series
// created by olav, Mi  2 Dez 2009 21:52:52 CET

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <ctime>
#include <cstring>
#include <sstream>
// #include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>

#include "../../Simulationen/olav.h"

#define REPORTS 25

#define AUTOREGRESSION_ORDER 3

// DemianTest with 10ms sampling
// #define NUM_NEURONS 100
// #define NUM_SAMPLES 126732
// #define DATA_BINS 15
// #define INPUTFILE "output/xresponse_10ms_15bins.dat"
// #define OUTPUTFILE "granger/grangercausality_os_10ms_15bins.mx"

// DemianTest with 20ms sampling
#define NUM_NEURONS 100
#define NUM_SAMPLES 89800
// binnning information is only relevant for comparison of to GC evaluations:
// #define DATA_BINS 15
#define INPUTFILE "output/xresponse_15bins.dat"
#define OUTPUTFILE "granger/grangercausality_os_15bins_p3.mx"



using namespace std;

void write_result(double **array);
void write_multidim_result(double ***array);
void load_data(double **array);

int main(int argc, char *argv[])
{
  cout <<"------ granger:main ------ olav, Wed 02 Dec 2009 ------"<<endl;
  time_t start, end;
  time_t middle;
  // unsigned long totaltrials, completedtrials;

  time(&start);
  cout <<"start: "<<ctime(&start)<<flush;
  cout <<"running on host: "<<flush;
  system("hostname");
  cout <<"current directory: "<<flush;
  system("pwd");

  cout <<"input file: "<<INPUTFILE<<endl;
  cout <<"output file: "<<OUTPUTFILE<<endl;

  cout <<"allocating memory..."<<flush;
  double **xdata = new double*[NUM_NEURONS];
  double ***xresult = new double**[NUM_NEURONS];
  for(int i=0; i<NUM_NEURONS; i++)
    xdata[i] = new double[NUM_SAMPLES];
  for(int i=0; i<NUM_NEURONS; i++)
	{
    xresult[i] = new double*[NUM_NEURONS];
    for(int j=0; j<NUM_NEURONS; j++)
			xresult[i][j] = new double[AUTOREGRESSION_ORDER];
  }

	// allocate GSL workspace
	gsl_multifit_linear_workspace * GSLworkspace = \
		gsl_multifit_linear_alloc(NUM_SAMPLES-AUTOREGRESSION_ORDER,2*AUTOREGRESSION_ORDER+1);
	double residue;
	gsl_matrix* input = gsl_matrix_alloc(NUM_SAMPLES-AUTOREGRESSION_ORDER,2*AUTOREGRESSION_ORDER+1);
	// because the last column of this matrix never changes, we can set it here:
	gsl_matrix_set_all(input,1.0);
	gsl_matrix* cov = gsl_matrix_alloc(2*AUTOREGRESSION_ORDER+1,2*AUTOREGRESSION_ORDER+1);
	gsl_vector* output = gsl_vector_alloc(NUM_SAMPLES-AUTOREGRESSION_ORDER);
	gsl_vector* coeff = gsl_vector_alloc(2*AUTOREGRESSION_ORDER+1);

  cout <<" done."<<endl;

  cout <<"loading data..."<<flush;
  load_data(xdata);
  cout <<" done."<<endl;

  // main loop:
	cout <<"set-up: "<<NUM_NEURONS<<" neurons";
	cout <<", "<<NUM_NEURONS*(NUM_NEURONS-1)<<" trials";
	cout <<", regression order "<<AUTOREGRESSION_ORDER;
	cout <<", "<<NUM_SAMPLES<<" samples"<<endl;

  for(int ii=0; ii<NUM_NEURONS; ii++)
  {
		cout <<"node #"<<ii+1<<": "<<flush;

		// set up source data of connection
		for(long tt=0; tt<NUM_SAMPLES-AUTOREGRESSION_ORDER; tt++)
			for(long tt2=0; tt2<AUTOREGRESSION_ORDER; tt2++)
				gsl_matrix_set(input,tt,AUTOREGRESSION_ORDER+tt2,xdata[ii][tt-tt2]);
		
    for(int jj=0; jj<NUM_NEURONS; jj++)
    {
	  	status(jj,REPORTS,NUM_NEURONS);
      if (ii != jj)
      {
				// add target data of connection and perform auto-regression
				for(long tt=0; tt<NUM_SAMPLES-AUTOREGRESSION_ORDER; tt++)
				{
					gsl_vector_set(output,tt,xdata[jj][tt+1]);
					
					for(long tt2=0; tt2<AUTOREGRESSION_ORDER; tt2++)
						gsl_matrix_set(input,tt,tt2,xdata[jj][tt-tt2]);
				}
				
				gsl_multifit_linear(input,output,coeff,cov,&residue,GSLworkspace);
				// save those coefficients that depend on the source
				for (int kk=0; kk<AUTOREGRESSION_ORDER; kk++)
					xresult[ii][jj][kk] = gsl_vector_get(coeff,AUTOREGRESSION_ORDER+kk);
      }
      else for (int kk=0; kk<AUTOREGRESSION_ORDER; kk++)
				xresult[ii][jj][kk] = 0.0;
    }
		time(&middle);
		cout <<" (elapsed "<<sec2string(difftime(middle,start))<<",";
		cout <<" ETA "<<ETAstring(ii+1,NUM_NEURONS,difftime(middle,start))<<")"<<endl;
  }

  time(&end);
  cout <<"end: "<<ctime(&end)<<flush;
  cout <<"runtime: "<<sec2string(difftime(end,start))<<endl;

  write_multidim_result(xresult);

  // free memory
	gsl_multifit_linear_free(GSLworkspace);
	gsl_matrix_free(input);
	gsl_matrix_free(cov);
	gsl_vector_free(output);
	gsl_vector_free(coeff);
	for(int i=0; i<NUM_NEURONS; i++) delete[] xdata[i];
  delete[] xdata;
  for(int i=0; i<NUM_NEURONS; i++) delete[] xresult[i];
  delete[] xresult;

  return 0;
}

void write_result(double **array)
{
  ofstream fileout1(OUTPUTFILE);

  bool already;
  fileout1 <<"{";
  for(int j=0; j<NUM_NEURONS; j++)
  {
  	if(j>0) fileout1<<",";
  	fileout1 <<"{";
  	already = false;
    for(unsigned long i=0; i<NUM_NEURONS; i++)
    {
      if (already) fileout1<<",";
    	fileout1 <<(double)array[j][i];
    	already = true;
    }
    fileout1 <<"}"<<endl;
  }
  fileout1 <<"}"<<endl;

  cout <<"Granger causality matrix saved."<<endl;
}
void write_multidim_result(double ***array)
{
  ofstream fileout1(OUTPUTFILE);

  fileout1 <<"{";
  for(int j=0; j<NUM_NEURONS; j++)
  {
  	if(j>0) fileout1<<",";
  	fileout1 <<"{";
    for(unsigned long i=0; i<NUM_NEURONS; i++)
    {
      if (i>0) fileout1<<",";
			fileout1 <<"{";
	    for(int k=0; k<AUTOREGRESSION_ORDER; k++)
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
}

void load_data(double **array)
{
  ifstream binaryfile(INPUTFILE, ios::binary);
	char* temparray = new char[NUM_SAMPLES];

  if (binaryfile == NULL)
  {
  	cout <<endl<<endl<<"error: cannot find input file! exiting."<<endl;
  	exit(1);
  }

  for(int j=0; j<NUM_NEURONS; j++)
  {
    binaryfile.read(temparray, NUM_SAMPLES);

		if (binaryfile.eof())
		{
			cout <<"Error: end of input file!"<<endl;
			exit(1);
		}

    for(long k=0; k<NUM_SAMPLES; k++)
			array[j][k] = double(temparray[k]);
  }

	if (binaryfile.peek() != EOF)
		cout <<"Warning: input file not completely read, parameters may be wrong."<<endl;
}
