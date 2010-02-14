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
#undef SHOW_DETAILED_PROGRESS

// currently this only works for a word length of 1
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
#define DATA_BINS 5
#define INPUTFILE "test/xresponse_5bins.dat"
#define OUTPUTFILE "test/transferentropy_os_15bins_global_test.mx"

using namespace std;

double TransferEntropy(char *arrayI, char *arrayJ, char *global, unsigned long** F_Ipast_Gpast, unsigned long*** F_Inow_Ipast_Gpast, unsigned long*** F_Ipast_Jpast_Gpast, unsigned long**** F_Inow_Ipast_Jpast_Gpast);
void write_result(double **array);
void load_data(char **array);
// bool match_backwards(char* x, unsigned long x_offset, char* y, unsigned long y_offset, unsigned int length);
// bool next_char(char* vector);
// bool next_char(char* vector, unsigned const int length);
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

	// allocate matrices for transition frequency of TE calculation
	unsigned long** F_Ipast_Gpast = new unsigned long*[DATA_BINS];
	unsigned long*** F_Inow_Ipast_Gpast = new unsigned long**[DATA_BINS];
	unsigned long*** F_Ipast_Jpast_Gpast = new unsigned long**[DATA_BINS];
	unsigned long**** F_Inow_Ipast_Jpast_Gpast = new unsigned long***[DATA_BINS];
	for (char x=0; x<DATA_BINS; x++)
	{
		F_Ipast_Gpast[x] = new unsigned long[DATA_BINS];
		F_Inow_Ipast_Gpast[x] = new unsigned long*[DATA_BINS];
		F_Ipast_Jpast_Gpast[x] = new unsigned long*[DATA_BINS];
		F_Inow_Ipast_Jpast_Gpast[x] = new unsigned long**[DATA_BINS];
		for (char x2=0; x2<DATA_BINS; x2++)
		{
			F_Inow_Ipast_Gpast[x][x2] = new unsigned long[DATA_BINS];
			F_Ipast_Jpast_Gpast[x][x2] = new unsigned long[DATA_BINS];
			F_Inow_Ipast_Jpast_Gpast[x][x2] = new unsigned long*[DATA_BINS];
			for (char x3=0; x3<DATA_BINS; x3++)
				F_Inow_Ipast_Jpast_Gpast[x][x2][x3] = new unsigned long[DATA_BINS];
		}
	}	

  cout <<" done."<<endl;

  cout <<"loading data..."<<flush;
  load_data(xdata);
  cout <<" done."<<endl;

  cout <<"generating global signal..."<<flush;
	char* xglobal = new char[NUM_SAMPLES];
  generate_global(xdata,xglobal);
	// randomize_signal(xglobal,5);
  cout <<" done."<<endl;
	
  // main loop:
  totaltrials = NUM_NEURONS*(NUM_NEURONS-1);
  cout <<"set-up: "<<NUM_NEURONS<<" neurons, "<<NUM_SAMPLES<<" samples, "<<DATA_BINS<<" bins"<<endl;
	cout <<"assumed length of Markov chain: "<<WORD_LENGTH<<endl;
  completedtrials = 0;
	// unsigned long long terms_sum = 0;
	// unsigned long long terms_zero = 0;
		
  for(int ii=0; ii<NUM_NEURONS; ii++)
  {
#ifndef SHOW_DETAILED_PROGRESS
  	cout <<"#"<<ii+1<<": "<<flush;
#endif
    for(int jj=0; jj<NUM_NEURONS; jj++)
    {
#ifndef SHOW_DETAILED_PROGRESS
  		// status(jj,REPORTS,NUM_NEURONS);
#endif
      if (ii != jj)
      {
#ifdef SHOW_DETAILED_PROGRESS
      	cout <<"#"<<ii+1<<" -> #"<<jj+1<<": "<<flush;
#endif
      	xresult[ii][jj] = TransferEntropy(xdata[jj], xdata[ii], xglobal, F_Ipast_Gpast, F_Inow_Ipast_Gpast, F_Ipast_Jpast_Gpast, F_Inow_Ipast_Jpast_Gpast);
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
		cout <<"done. (elapsed "<<sec2string(difftime(middle,start))<<",";
		cout <<" ETA "<<ETAstring(ii+1,NUM_NEURONS,difftime(middle,start))<<")"<<endl;
#endif
  }
#ifndef SHOW_DETAILED_PROGRESS
  cout <<" done."<<endl;
#endif

  time(&end);
  cout <<"end: "<<ctime(&end)<<flush;
  cout <<"runtime: "<<sec2string(difftime(end,start))<<endl;

	// cout <<"TE terms: "<<terms_sum<<", of those zero: "<<terms_zero<<" ("<<int(double(terms_zero)*100/terms_sum)<<"%)"<<endl;

  write_result(xresult);

  for(int i=0; i<NUM_NEURONS; i++) delete[] xdata[i];
  delete[] xdata;
  for(int i=0; i<NUM_NEURONS; i++) delete[] xresult[i];
  delete[] xresult;

  return 0;
}

double TransferEntropy(char *arrayI, char *arrayJ, char *global, unsigned long** F_Ipast_Gpast, unsigned long*** F_Inow_Ipast_Gpast, unsigned long*** F_Ipast_Jpast_Gpast, unsigned long**** F_Inow_Ipast_Jpast_Gpast)
{
	/* see for reference:
	     Gourevitch und Eggermont. Evaluating Information Transfer Between Auditory
	     Cortical Neurons. Journal of Neurophysiology (2007) vol. 97 (3) pp. 2533 */
  // double result = 0.0;

	// We are looking at the information flow of arrayJ -> arrayI
	
	// clear memory
	for (char x=0; x<DATA_BINS; x++)
	{
		memset(F_Ipast_Gpast[x], 0, DATA_BINS*sizeof(unsigned long));
		for (char x2=0; x2<DATA_BINS; x2++)
		{
			memset(F_Inow_Ipast_Gpast[x][x2], 0, DATA_BINS*sizeof(unsigned long));
			memset(F_Ipast_Jpast_Gpast[x][x2], 0, DATA_BINS*sizeof(unsigned long));
			for (char x3=0; x3<DATA_BINS; x3++)
				memset(F_Inow_Ipast_Jpast_Gpast[x][x2][x3], 0, DATA_BINS*sizeof(unsigned long));
		}
	}
	
  // extract probabilities (actually number of occurrence)
	for (unsigned long t=WORD_LENGTH; t<NUM_SAMPLES; t++)
  {
		F_Ipast_Gpast[arrayI[t-1]][global[t-1]]++;
		F_Inow_Ipast_Gpast[arrayI[t]][arrayI[t-1]][global[t-1]]++;
		F_Ipast_Jpast_Gpast[arrayI[t-1]][arrayJ[t-1]][global[t-1]]++;
		F_Inow_Ipast_Jpast_Gpast[arrayI[t]][arrayI[t-1]][arrayJ[t-1]][global[t-1]]++;
#ifdef SHOW_DETAILED_PROGRESS
		status(t, REPORTS, NUM_SAMPLES-WORD_LENGTH);
#endif
	}
	
	// index convention:
	// k - Ipast
	// l - Jpast
	// m - Inow
	// g - Gpast
	
	// calculate transfer entropy
	// for (char k=0; k<DATA_BINS; k++)
	// 	for (char g=0; g<DATA_BINS; g++)
	// 		if (F_Ipast_Gpast[k][g]!=0) for (char l=0; l<DATA_BINS; l++)
	// 			if (F_Ipast_Jpast_Gpast[k][l][g]!=0) for (char m=0; m<DATA_BINS; m++)
	// 					if ((F_Inow_Ipast_Gpast[m][k][g]!=0) && (F_Inow_Ipast_Jpast_Gpast[m][k][l][g] != 0))
	// 						result += double(F_Inow_Ipast_Jpast_Gpast[m][k][l][g])/double(NUM_SAMPLES-WORD_LENGTH) * log(double(F_Inow_Ipast_Jpast_Gpast[m][k][l][g])*double(F_Ipast_Gpast[k][g])/(double(F_Ipast_Jpast_Gpast[k][l][g])*double(F_Inow_Ipast_Gpast[m][k][g])));
	double Hxx = 0.0;
	for (char k=0; k<DATA_BINS; k++)
		for (char m=0; m<DATA_BINS; m++)
			for (char g=0; g<DATA_BINS; g++)
				if (F_Inow_Ipast_Gpast[m][k][g] > 0)
					Hxx -= double(F_Inow_Ipast_Gpast[m][k][g])/(NUM_SAMPLES-WORD_LENGTH) * log(double(F_Inow_Ipast_Gpast[m][k][g])/double(F_Ipast_Gpast[k][g]));
	Hxx /= log(2);
	
	double Hxxy = 0.0;
	for (char k=0; k<DATA_BINS; k++)
		for (char l=0; l<DATA_BINS; l++)
			for (char m=0; m<DATA_BINS; m++)
				for (char g=0; g<DATA_BINS; g++)
					if (F_Inow_Ipast_Jpast_Gpast[m][k][l][g] > 0)
						Hxxy -= double(F_Inow_Ipast_Jpast_Gpast[m][k][l][g])/(NUM_SAMPLES-WORD_LENGTH) * log(double(F_Inow_Ipast_Jpast_Gpast[m][k][l][g])/double(F_Ipast_Jpast_Gpast[k][l][g]));
	Hxxy /= log(2);

  return (Hxx - Hxxy);
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

/* INLINE bool match_backwards(char* x, unsigned long x_offset, char* y, unsigned long y_offset, unsigned int length)
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
} */

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

/* unsigned long index(char v1)
{
	return (unsigned long)(v1);
}

unsigned long index(char v1, char v2)
{
	return (unsigned long)(v1)*DATA_BINS + v2;
}

unsigned long index(char v1, char v2, char v3)
{
	return (unsigned long)(v1)*DATA_BINS*DATA_BINS + (unsigned long)(v2)*DATA_BINS + v3;
} */
