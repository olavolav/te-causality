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

#define AUTOREGRESSION_ORDER 1
#define RESULT_DIMENSION AUTOREGRESSION_ORDER

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
#define INPUTFILE "test/xresponse_15bins.dat"
#define OUTPUTFILE "test/grangercausality_os_15bins-allatoncetest.mx"



using namespace std;

void write_result(double **array);
void write_multidim_result(double ***array);
void load_data(double **array);

int main(int argc, char *argv[])
{
  cout <<"------ granger:allatonce ------ olav, Sun 24 Jan 2010 ------"<<endl;
  time_t start, end;
  time_t middle;

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
  for(int i=0; i<NUM_NEURONS; i++)
    xdata[i] = new double[NUM_SAMPLES];

#if RESULT_DIMENSION > 1
  double ***xresult = new double**[NUM_NEURONS];
  for(int i=0; i<NUM_NEURONS; i++)
	{
    xresult[i] = new double*[NUM_NEURONS];
    for(int j=0; j<NUM_NEURONS; j++)
			xresult[i][j] = new double[AUTOREGRESSION_ORDER+1];
  }
#else
  double **xresult = new double*[NUM_NEURONS];
  for(int i=0; i<NUM_NEURONS; i++)
    xresult[i] = new double[NUM_NEURONS];
#endif

	// allocate GSL workspaces
	gsl_multifit_linear_workspace* GSLworkspaceBoth = \
		gsl_multifit_linear_alloc(NUM_SAMPLES-AUTOREGRESSION_ORDER,NUM_NEURONS*AUTOREGRESSION_ORDER+1);
	double residue;
	// the coulmns in input are: (nr. of columns)
	// target (AUTOREGRESSION_ORDER), source (AUTOREGRESSION_ORDER), constant term (1)
	gsl_matrix* inputBoth = gsl_matrix_alloc(NUM_SAMPLES-AUTOREGRESSION_ORDER,NUM_NEURONS*AUTOREGRESSION_ORDER+1);
	// because the last column of this matrix never changes, we can set it here:
	gsl_matrix_set_all(inputBoth,1.0);
	gsl_matrix* covBoth = gsl_matrix_alloc(NUM_NEURONS*AUTOREGRESSION_ORDER+1,NUM_NEURONS*AUTOREGRESSION_ORDER+1);
	// output is the only GSL consturuct that we use for both cases ("both" and "single")
	gsl_vector* output = gsl_vector_alloc(NUM_SAMPLES-AUTOREGRESSION_ORDER);
	gsl_vector* coeffBoth = gsl_vector_alloc(NUM_NEURONS*AUTOREGRESSION_ORDER+1);

	// gsl_multifit_linear_workspace * GSLworkspaceSingle =	gsl_multifit_linear_alloc(NUM_SAMPLES-AUTOREGRESSION_ORDER,1*AUTOREGRESSION_ORDER+1);
	// gsl_matrix* inputSingle = gsl_matrix_alloc(NUM_SAMPLES-AUTOREGRESSION_ORDER,1*AUTOREGRESSION_ORDER+1);
	// gsl_matrix_set_all(inputSingle,1.0);
	// gsl_matrix* covSingle = gsl_matrix_alloc(1*AUTOREGRESSION_ORDER+1,1*AUTOREGRESSION_ORDER+1);
	// gsl_vector* coeffSingle = gsl_vector_alloc(1*AUTOREGRESSION_ORDER+1);

  cout <<" done."<<endl;

  cout <<"loading data..."<<flush;
  load_data(xdata);
  cout <<" done."<<endl;

	cout <<"set-up: "<<NUM_NEURONS<<" neurons";
	// cout <<", "<<NUM_NEURONS*(NUM_NEURONS-1)<<" trials";
	cout <<", regression order "<<AUTOREGRESSION_ORDER;
	cout <<", "<<NUM_SAMPLES<<" samples"<<endl;
	
	cout <<"initializing workspace (input side)..."<<flush;
	// we can set it here since it is the same for all neurons
	for (int tt=AUTOREGRESSION_ORDER; tt<NUM_SAMPLES; tt++)
		for (int jj=0; jj<NUM_NEURONS; jj++)
			for (int nn=0; nn<AUTOREGRESSION_ORDER; nn++)
				gsl_matrix_set(inputBoth,tt-AUTOREGRESSION_ORDER,AUTOREGRESSION_ORDER*jj+nn,xdata[jj][tt-nn-1]);
 	cout <<" done."<<endl;

	// double test_residue;
	// double covTraceBoth, covTraceSingle;
  // main loop:
  for(int ii=0; ii<NUM_NEURONS; ii++)
  {
  	// status(ii,REPORTS,NUM_NEURONS);
		cout <<"#"<<ii+1<<"... "<<flush;

		// add target data of connection
		for(long tt=AUTOREGRESSION_ORDER; tt<NUM_SAMPLES; tt++)
			gsl_vector_set(output,tt-AUTOREGRESSION_ORDER,xdata[ii][tt]);
		
		// perform auto-regression
		gsl_multifit_linear(inputBoth,output,coeffBoth,covBoth,&residue,GSLworkspaceBoth);
		
		// test_residue = residue;
		// gsl_multifit_linear(inputSingle,output,coeffSingle,covSingle,&residue,GSLworkspaceSingle);

		// test:
		// xresult[ii][jj] = residue/test_residue;
		
		// save those coefficients that depend on the source
		// covTraceBoth = covTraceSingle = 0.0;
		// for (int kk=0; kk<AUTOREGRESSION_ORDER; kk++)
		// {
		// 	xresult[ii][jj][kk] = gsl_vector_get(coeffBoth,AUTOREGRESSION_ORDER+kk);
			// covTraceBoth += gsl_matrix_get(covBoth,AUTOREGRESSION_ORDER+kk,AUTOREGRESSION_ORDER+kk);
			// covTraceSingle += gsl_matrix_get(covSingle,AUTOREGRESSION_ORDER+kk,AUTOREGRESSION_ORDER+kk);
		// }
		// xresult[ii][jj][AUTOREGRESSION_ORDER] = residue;
	  for(int jj=0; jj<NUM_NEURONS; jj++)
			xresult[jj][ii] = gsl_vector_get(coeffBoth,ii*AUTOREGRESSION_ORDER);
		
		// save redultion in variance as xresult
		// xresult[ii][jj] = log(sqrt(gsl_matrix_get(covSingle,0,0)+gsl_matrix_get(covSingle,1,1))/sqrt(gsl_matrix_get(covBoth,0,0)+gsl_matrix_get(covBoth,1,1)));
		time(&middle);
		cout <<" (elapsed "<<sec2string(difftime(middle,start))<<",";
		cout <<" ETA "<<ETAstring(ii+1,NUM_NEURONS,difftime(middle,start))<<")"<<endl;
  }

  time(&end);
  cout <<"end: "<<ctime(&end)<<flush;
  cout <<"runtime: "<<sec2string(difftime(end,start))<<endl;

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
	for(int i=0; i<NUM_NEURONS; i++) delete[] xdata[i];
  delete[] xdata;
  for(int i=0; i<NUM_NEURONS; i++) delete[] xresult[i];
  delete[] xresult;

  return 0;
}

void write_result(double **array)
{
  ofstream fileout1(OUTPUTFILE);

	fileout1.precision(6);
	fileout1 <<fixed;
  fileout1 <<"{";
  for(int j=0; j<NUM_NEURONS; j++)
  {
  	if(j>0) fileout1<<",";
  	fileout1 <<"{";
    for(unsigned long i=0; i<NUM_NEURONS; i++)
    {
      if (i>0) fileout1<<",";
    	fileout1 <<(double)array[j][i];
    }
    fileout1 <<"}"<<endl;
  }
  fileout1 <<"}"<<endl;

  cout <<"Granger causality matrix saved."<<endl;
}

void write_multidim_result(double ***array)
{
  ofstream fileout1(OUTPUTFILE);

	fileout1.precision(9);
	fileout1 <<fixed;
  fileout1 <<"{";
  for(unsigned long j=0; j<NUM_NEURONS; j++)
  {
  	if(j>0) fileout1<<",";
  	fileout1 <<"{";
    for(unsigned long i=0; i<NUM_NEURONS; i++)
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
}

void load_data(double **array)
{
  ifstream binaryfile(INPUTFILE, ios::binary);
	char* temparray = new char[NUM_SAMPLES];

  if (binaryfile == NULL)
  {
  	cout <<endl<<"error: cannot find input file! exiting."<<endl;
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
