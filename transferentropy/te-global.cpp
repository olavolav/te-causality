// calculate the transfer entropy between a numer of time series (including global signal)
// created by olav, Mon 18 Jan 2010

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
// #define NUM_SAMPLES 126732
// DemianTest with 20ms sampling
#define NUM_NEURONS 100
#define NUM_SAMPLES 89800
#define DATA_BINS 15
#define INPUTFILE "test/xresponse_15bins.dat"
#define OUTPUTFILE "test/transferentropy_os_15bins_global.mx"

using namespace std;

double TEterm(char *array1, char *array2, char *global, char k, char* l, char* m, unsigned long long* terms_sum, unsigned long long* terms_zero);
double TransferEntropy(char *array1, char *array2, char *global, unsigned long long* terms_sum, unsigned long long* terms_zero);
void write_result(double **array);
void load_data(char **array);
bool match_backwards(char* x, unsigned long x_offset, char* y, unsigned long y_offset, unsigned int length);
bool next_char(char* vector);
bool next_char(char* vector, unsigned const int length);
void generate_global(char** raw, char* global);
void randomize_signal(char* raw, char deviation_in_bins);

int main(int argc, char *argv[])
{
  cout <<"------ transferentropy:global ------ olav, Mon 18 Jan 2010 ------"<<endl;
  time_t start, end, middle;
  unsigned long totaltrials, completedtrials;

  time(&start);
  cout <<"start: "<<ctime(&start)<<flush;
  cout <<"running on host: "<<flush;
  system("hostname");
  cout <<"current directory: "<<flush;
  system("pwd");

  cout <<"input file: "<<INPUTFILE<<endl;
  cout <<"output file: "<<OUTPUTFILE<<endl;

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

  cout <<"generating global signal..."<<flush;
	char* xglobal = new char[NUM_SAMPLES];
  generate_global(xdata,xglobal);
  cout <<" done."<<endl;
  // cout <<"randomizing global signal..."<<flush;
  // randomize_signal(xglobal,10);
  // cout <<" done."<<endl;

	// for (int t=1; t<25; t++)
	// {
	// 	cout <<"t="<<t<<": local_3 "<<int(xdata[3][t])<<", global "<<int(xglobal[t]);
	// 	cout <<", match: "<<int(match_backwards(xdata[3],t,xglobal,t,WORD_LENGTH));
	// 	cout <<endl;
	// }
	// exit (1);
	
  // main loop:
  totaltrials = NUM_NEURONS*(NUM_NEURONS-1);
  cout <<"set-up: "<<NUM_NEURONS<<" neurons, "<<NUM_SAMPLES<<" samples, "<<DATA_BINS<<" bins"<<endl;
	cout <<"assumed length of Markov chain: "<<WORD_LENGTH<<endl;
  completedtrials = 0;
	unsigned long long terms_sum = 0;
	unsigned long long terms_zero = 0;
	
	// char* l = new char[2*WORD_LENGTH];
	// memset(l, 0, 2*WORD_LENGTH*sizeof(char));
	// do
	// {
	// 	cout <<"newline: ";
	// 	for (int i=0; i<2*WORD_LENGTH; i++)
	// 		cout <<int(l[i])<<" ";
	// 	cout <<endl;
	// } while (next_char(l,2*WORD_LENGTH));
	
  for(int ii=0; ii<NUM_NEURONS; ii++)
  {
#ifndef SHOW_DETAILED_PROGRESS
  	cout <<"#"<<ii+1<<": "<<flush;
#endif
    for(int jj=0; jj<NUM_NEURONS; jj++)
    {
#ifndef SHOW_DETAILED_PROGRESS
  		status(jj,REPORTS,NUM_NEURONS);
#endif
      if (ii != jj)
      {
#ifdef SHOW_DETAILED_PROGRESS
      	cout <<"#"<<ii+1<<" -> #"<<jj+1<<": "<<flush;
#endif
      	xresult[ii][jj] = TransferEntropy(xdata[ii], xdata[jj], xglobal, &terms_sum, &terms_zero);
				completedtrials++;
#ifdef SHOW_DETAILED_PROGRESS
				time(&middle);
				cout <<" (elapsed "<<sec2string(difftime(middle,start))<<",";
				cout <<" ETA "<<ETAstring(completedtrials,totaltrials,difftime(middle,start))<<")"<<endl;
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
#ifndef SHOW_DETAILED_PROGRESS
		time(&middle);
		cout <<" (elapsed "<<sec2string(difftime(middle,start))<<",";
		cout <<" ETA "<<ETAstring(ii+1,NUM_NEURONS,difftime(middle,start))<<")"<<endl;
#endif
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

INLINE double TEterm(char *array1, char *array2, char *global, char k, char* l, char* m, unsigned long long* terms_sum, unsigned long long* terms_zero)
{
  unsigned long countA, countB, countC, countD;
  countA = countB = countC = countD = 0;
  double result = 0.0;

	// relation to Schreiber et al. (J->I):
	//   k corresponds to i_{n+1}
	//   l corresponds to i_{n}^k
	//   m corresponds to j_{n}^l
	// therefore:
	//   countD ~ P(i_{n}^k)
	//   countC ~ P(i_{n+1}, i_{n}^k)
	//   countB ~ P(j_{n}^l, i_{n}^k)
	//   countA ~ P(i_{n+1}, j_{n}^l, i_{n}^k)
	
  for (unsigned long tt=WORD_LENGTH; tt<NUM_SAMPLES; tt++)
  {
		if (match_backwards(global,tt-1,l,2*WORD_LENGTH-1,WORD_LENGTH) && \
			match_backwards(array2,tt-1,l,WORD_LENGTH-1,WORD_LENGTH))
		{
  		countD++;
  		// old version:
			// if(array2[tt] == k) countC++;
  		// if(match_backwards(array1,tt-1,m,WORD_LENGTH-1,WORD_LENGTH))
  		// {
  		// 	countB++;
  		// 	if(array2[tt] == k) countA++;
  		// }

			// new version, marginally faster:
			if(array2[tt] == k)
			{
				countC++;
		  	if(match_backwards(array1,tt-1,m,WORD_LENGTH-1,WORD_LENGTH)) countA++;
			}
			else if(match_backwards(array1,tt-1,m,WORD_LENGTH-1,WORD_LENGTH)) countB++;
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


double TransferEntropy(char *array1, char *array2, char *global, unsigned long long* terms_sum, unsigned long long* terms_zero)
{
	/* see for reference:
	     Gourevitch und Eggermont. Evaluating Information Transfer Between Auditory
	     Cortical Neurons. Journal of Neurophysiology (2007) vol. 97 (3) pp. 2533 */
  double result = 0.0;

	// We are looking at the information flow of array1 -> array2
	
	// as opposed to te-main, l now includes the global signal:
	char* l = new char[2*WORD_LENGTH];
	memset(l, 0, 2*WORD_LENGTH*sizeof(char));
	char* m = new char[WORD_LENGTH];
	
#ifdef SHOW_DETAILED_PROGRESS
	int const max_index = int(pow(double(DATA_BINS),3*WORD_LENGTH));
	unsigned long running_index = 0;
#endif
	do
	{
		memset(m, 0, WORD_LENGTH*sizeof(char));
		do
    {
#ifdef SHOW_DETAILED_PROGRESS
			status(running_index, REPORTS, max_index);
			running_index++;
#endif
      for (char k=0; k<DATA_BINS; k++)
				// cout <<"debug: k=("<<int(k)<<"), l=("<<int(l[0])<<","<<int(l[1])<<"), m=("<<int(m[0])<<")"<<endl;
      	result += TEterm(array1, array2, global, k, l, m, terms_sum, terms_zero);
		}
		while (next_char(m));
	}
	while (next_char(l,2*WORD_LENGTH));

  return result;
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
	for (unsigned int i=0; i<length; i++)
		if (x[x_offset-i] != y[y_offset-i])
			return false;
	
	return true;
}

INLINE bool next_char(char* vector)
{
	return next_char(vector, WORD_LENGTH);
}

bool next_char(char* vector, unsigned const int length)
{
	bool not_at_end = false;
	
	for (unsigned int i=0; i<length; i++)
		if (vector[i] != DATA_BINS-1)
		{
			not_at_end = true;
			break;
		}
	
	vector[0]++;
	for (unsigned int i=0; i<length; i++)
	{
		if (vector[i] >= DATA_BINS)
		{
			vector[i] -= DATA_BINS;
			// vector[i] = 0;
			if (i+1 < length) vector[i+1]++;
		}
		else break;
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
		global[t] = char(round(avg/NUM_NEURONS));
		// test:
		// global[t] = 3;
	}
}

void randomize_signal(char* raw, char deviation_in_bins)
{
	int xtemp;
	
	for (unsigned long t=0; t<NUM_SAMPLES; t++)
	{
		xtemp = raw[t] + char((double(rand())/RAND_MAX)*2*deviation_in_bins-deviation_in_bins);
		if (xtemp<0) xtemp = 0;
		else if (xtemp>=DATA_BINS) xtemp = DATA_BINS-1;
		raw[t] = char(xtemp);
	}
}
