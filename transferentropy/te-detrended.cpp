// calculate the transfer entropy between a numer of time series (detrended version)
// created by olav, Do 14 Jan 2010 17:33:11 CET

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <ctime>
#include <cstring>
#include <sstream>

#include "../../Simulationen/olav.h"

#ifndef INLINE
#define INLINE extern inline
#endif

#define REPORTS 25
#define SHOW_DETAILED_PROGRESS

#define WORD_LENGTH 1

// #define NUM_NEURONS 50
// samples: 33331, 40177, 51025, 59777
// #define NUM_SAMPLES 59777
// #define DATA_BINS 30
// actual time series
// #define NUM_NEURONS 100
// #define NUM_SAMPLES 59951
// #define DATA_BINS 30
// actual time series, real data
// #define NUM_NEURONS 134
// #define NUM_SAMPLES 63536
// #define DATA_BINS 30

// maximum amplitudes only
// #define NUM_NEURONS 100
// #define NUM_SAMPLES 165
// #define DATA_BINS 15
// maximum amplitudes only, real data
// #define NUM_NEURONS 134
// #define NUM_SAMPLES 277
// #define DATA_BINS 20

// DemianTest with 10ms sampling
// #define NUM_NEURONS 100
// #define NUM_SAMPLES 126732
// #define DATA_BINS 15
// #define INPUTFILE "output/xresponse_10ms_15bins.dat"
// #define OUTPUTFILE "output/transferentropy_os_10ms_15bins.mx"
// DemianTest with 20ms sampling
#define NUM_NEURONS 100
#define NUM_SAMPLES 89800
#define DATA_BINS 30
#define DATA_BINS_ONLY_FOR_PADDING 15
#define INPUTFILE "output/xresponse_15bins.dat"
#define OUTPUTFILE "te/transferentropy_os_15bins_detrended-test_to30bins.mx"

using namespace std;

double TEterm(char *array1, char *array2, char k, char* l, char* m, unsigned long long* counter);
double TransferEntropy(char *array1, char *array2, unsigned long long* terms_sum, unsigned long long* terms_zero);
void write_result(double **array);
void load_data(char **array);
bool match_backwards(char* x, unsigned long x_offset, char* y, unsigned long y_offset, unsigned int length);
bool next_char(char* vector);
void generate_global(char** raw, char* global);
void detrend_signal(char** data, char* global);


int main(int argc, char *argv[])
{
  cout <<"------ transferentropy:detrended ------ olav, Thu 14 Jan 2010 ------"<<endl;
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

  cout <<"input file: "<<INPUTFILE<<endl;
  cout <<"output file: "<<OUTPUTFILE<<endl;

	/* cout <<"TEST:"<<endl;
	char* xtest = new char[WORD_LENGTH];
	memset(xtest, 0, WORD_LENGTH*sizeof(char));
	for(int i=0; i<15*15+2; i++)
	{
		for (int j=0; j<WORD_LENGTH; j++) cout <<int(xtest[j])<<" ";
		cout <<" -> "<<int(next_char(xtest))<<endl;
	} */

  cout <<"allocating memory..."<<flush;
  char **xdata = new char*[NUM_NEURONS];
  double **xresult = new double*[NUM_NEURONS];
  for(int i=0; i<NUM_NEURONS; i++)
  {
    xdata[i] = new char[NUM_SAMPLES];
    memset(xdata[i], 0, NUM_SAMPLES*sizeof(char));
    xresult[i] = new double[NUM_NEURONS];
    memset(xresult[i], 0, NUM_NEURONS*sizeof(double));
  }
  cout <<" done."<<endl;

  cout <<"loading data..."<<flush;
  load_data(xdata);
  cout <<" done."<<endl;

  cout <<"generating global signal and detrending local signals..."<<flush;
	char* xglobal = new char[NUM_SAMPLES];
  generate_global(xdata,xglobal);
	detrend_signal(xdata,xglobal);
  cout <<" done."<<endl;

  // main loop:
  totaltrials = NUM_NEURONS*(NUM_NEURONS-1);
  cout <<"set-up: "<<NUM_NEURONS<<" neurons, "<<totaltrials<<" trials"<<endl;
	cout <<"assumed length of Markov chain: "<<WORD_LENGTH<<endl;
  completedtrials = 0;
	unsigned long long terms_sum = 0;
	unsigned long long terms_zero = 0;
	
#ifndef SHOW_DETAILED_PROGRESS
  	cout <<"running "<<flush;
#endif

  for(int ii=0; ii<NUM_NEURONS; ii++)
  {
#ifndef SHOW_DETAILED_PROGRESS
  	status(ii,REPORTS,NUM_NEURONS);
#endif
    for(int jj=0; jj<NUM_NEURONS; jj++)
    {
      if (ii != jj)
      {
#ifdef SHOW_DETAILED_PROGRESS
      	cout <<"#"<<ii+1<<" -> #"<<jj+1<<": "<<flush;
#endif
      	xresult[ii][jj] = TransferEntropy(xdata[ii], xdata[jj], &terms_sum, &terms_zero);
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

	cout <<"TE terms: "<<terms_sum<<", of those zero: "<<terms_zero<<" ("<<int(double(terms_zero)*100/terms_sum)<<"%)"<<endl;

  write_result(xresult);

  for(int i=0; i<NUM_NEURONS; i++) delete[] xdata[i];
  delete[] xdata;
  for(int i=0; i<NUM_NEURONS; i++) delete[] xresult[i];
  delete[] xresult;

  return 0;
}

INLINE double TEterm(char *array1, char *array2, char k, char* l, char* m, unsigned long long* terms_sum, unsigned long long* terms_zero)
{
  unsigned long countA, countB, countC, countD;
  countA = countB = countC = countD = 0;
  double result = 0.0;
	
	const char first_indication = l[WORD_LENGTH-1];

  for (unsigned long tt=WORD_LENGTH; tt<NUM_SAMPLES; tt++)
  {
		if (array2[tt-1] == first_indication) // optimization? (a few percent)
		{
#if WORD_LENGTH>1
			// if (match_backwards(array2,tt-1,l,WORD_LENGTH-1,WORD_LENGTH-1))
			if (match_backwards(array2,tt-1-1,l,WORD_LENGTH-1-1,WORD_LENGTH-1))
	  	{
#endif
	  		countD++;
	  		if(array2[tt] == k) countC++;
	  		if(match_backwards(array1,tt-1,m,WORD_LENGTH-1,WORD_LENGTH))
	  		{
	  			countB++;
	  			if(array2[tt] == k) countA++;
	  		}
#if WORD_LENGTH>1
	  	}
#endif
		}
  }

  if (countA*countB*countC*countD != 0)
  {
  	result = double(countA)/NUM_SAMPLES * log(double(countA*countD)/(countB*countC));

  	// transform unit of information to bits
  	result /= log(2);
  }
  else (*terms_zero)++;

	(*terms_sum)++;
	return result;
}

double TransferEntropy(char *array1, char *array2, unsigned long long* terms_sum, unsigned long long* terms_zero)
{
	/* see for reference:
	     Gourevitch und Eggermont. Evaluating Information Transfer Between Auditory
	     Cortical Neurons. Journal of Neurophysiology (2007) vol. 97 (3) pp. 2533 */
  double result = 0.0;

	char* l = new char[WORD_LENGTH];
	memset(l, 0, WORD_LENGTH*sizeof(char));
	char* m = new char[WORD_LENGTH];
		
	int const max_index = pow(double(DATA_BINS),2*WORD_LENGTH);
	unsigned long running_index = 0;
	do
	{
		memset(m, 0, WORD_LENGTH*sizeof(char));
		do
    {
#ifdef SHOW_DETAILED_PROGRESS
			status(running_index, REPORTS, max_index);
#endif
      for (char k=0; k<DATA_BINS; k++)
      	result += TEterm(array1, array2, k, l, m, terms_sum, terms_zero);
			running_index++;
		}
		while (next_char(m));
	}
	while (next_char(l));

  return result;
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
    fileout1 <<"}";
  }
  fileout1 <<"}"<<endl;

  cout <<"transfer entropy matrix saved."<<endl;
}

void load_data(char **array)
{
  ifstream binaryfile(INPUTFILE, ios::binary);

  if (binaryfile == NULL)
  {
  	cout <<endl<<endl<<"error: cannot find input file! exiting."<<endl;
  	exit(1);
  }

  for(int j=0; j<NUM_NEURONS; j++)
  {
    // status(j, REPORTS, NUM_NEURONS);
    binaryfile.read(array[j], NUM_SAMPLES);
  }

	if (binaryfile.eof())
	{
		cout <<"Error: end of input file!"<<endl;
		exit(1);
	}
	if (binaryfile.peek() != EOF)
		cout <<"Warning: input file not completely read, parameters may be wrong."<<endl;
}

INLINE bool match_backwards(char* x, unsigned long x_offset, char* y, unsigned long y_offset, unsigned int length)
{
	bool result = true;
	for (unsigned int i=0; i<length; i++)
		if (x[x_offset-i] != y[y_offset-i])
		{
			result = false;
			break;
		}
	
	return result;
}

bool next_char(char* vector)
{
	bool not_at_end = true;
	
	for (int i=0; i<WORD_LENGTH; i++)
		if (vector[i] != DATA_BINS-1) not_at_end = false;
	
	vector[0]++;
	for (int i=0; i<WORD_LENGTH; i++)
	{
		if (vector[i] >= DATA_BINS)
		{
			vector[i] -= DATA_BINS;
			if (i+1 < WORD_LENGTH) vector[i+1]++;
			// else not_at_end = false;
		}
	}
	
	return not_at_end;
}

void generate_global(char** raw, char* global)
{
	double avg;
	for (unsigned long t=0; t<NUM_SAMPLES; t++)
	{
		avg = 0.0;
		for (unsigned long j=0; j<NUM_NEURONS; j++)
			avg += double(raw[j][t]);
		global[t] = char(avg/NUM_NEURONS);
	}
}

void detrend_signal(char** data, char* global)
{
	int xtemp;
	for (unsigned long j=0; j<NUM_NEURONS; j++)
		for (unsigned long t=0; t<NUM_SAMPLES; t++)
		{
			xtemp = data[j][t] - global[t] + DATA_BINS_ONLY_FOR_PADDING;
			if (xtemp<0) xtemp = 0;
			else if (xtemp>=DATA_BINS) xtemp = DATA_BINS-1;
			data[j][t] = char(xtemp);
		}
}
