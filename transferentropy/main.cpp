// calculate the transfer entropy between a numer of time series
// created by olav, Mi 10 Jun 2009 19:42:11 IDT

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

// this version only works with 1!
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
// #define NUM_SAMPLES 87501
#define NUM_SAMPLES 89800
#define DATA_BINS 15
#define INPUTFILE "test/xresponse_15bins.dat"
// #define INPUTFILE "../../Simulationen/Imaging/DemianTest/FPSoutput2/noisefree_uchar_20ms_15bins.dat"
#define OUTPUTFILE "../transferentropy-sim/output/adjA_iteration2.mx"

using namespace std;

// double TEterm(char *array1, char *array2);
double TransferEntropy(char *arrayI, char *arrayJ, unsigned long* F_Ipast, unsigned long** F_Inow_Ipast, unsigned long** F_Ipast_Jpast, unsigned long*** F_Inow_Ipast_Jpast);
void write_result(double **array);
void load_data(char **array);
// bool match_backwards(char* x, unsigned long x_offset, char* y, unsigned long y_offset, unsigned int length);
// bool next_char(char* vector);
void generate_global(char** raw, char* global);

int main(int argc, char *argv[])
{
  cout <<"------ transferentropy:main ------ olav, Wed 10 Jun 2009 ------"<<endl;
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

	unsigned long* F_Ipast = new unsigned long[DATA_BINS];
	unsigned long** F_Inow_Ipast = new unsigned long*[DATA_BINS];
	unsigned long** F_Ipast_Jpast = new unsigned long*[DATA_BINS];
	unsigned long*** F_Inow_Ipast_Jpast = new unsigned long**[DATA_BINS];

	for (char x=0; x<DATA_BINS; x++)
	{
		F_Inow_Ipast[x] = new unsigned long[DATA_BINS];
		F_Ipast_Jpast[x] = new unsigned long[DATA_BINS];

		F_Inow_Ipast_Jpast[x] = new unsigned long*[DATA_BINS];
		for (char x2=0; x2<DATA_BINS; x2++)
		{
			F_Inow_Ipast_Jpast[x][x2] = new unsigned long[DATA_BINS];
		}
	}
  cout <<" done."<<endl;

  cout <<"loading data... "<<flush;
  load_data(xdata);
  cout <<" done."<<endl;

  // main loop:
  totaltrials = NUM_NEURONS*(NUM_NEURONS-1);
  cout <<"set-up: "<<NUM_NEURONS<<" neurons, "<<NUM_SAMPLES<<" samples, "<<DATA_BINS<<" bins"<<endl;
	cout <<"assumed length of Markov chain: "<<WORD_LENGTH<<endl;
  completedtrials = 0;
	// unsigned long long terms_sum = 0;
	// unsigned long long terms_zero = 0;
	
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
      	xresult[ii][jj] = TransferEntropy(xdata[ii], xdata[jj],F_Ipast, F_Inow_Ipast, F_Ipast_Jpast, F_Inow_Ipast_Jpast);
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

  write_result(xresult);

  for(int i=0; i<NUM_NEURONS; i++) delete[] xdata[i];
  delete[] xdata;
  for(int i=0; i<NUM_NEURONS; i++) delete[] xresult[i];
  delete[] xresult;

  return 0;
}

double TransferEntropy(char *arrayI, char *arrayJ, unsigned long* F_Ipast, unsigned long** F_Inow_Ipast, unsigned long** F_Ipast_Jpast, unsigned long*** F_Inow_Ipast_Jpast)
{
	/* see for reference:
	     Gourevitch und Eggermont. Evaluating Information Transfer Between Auditory
	     Cortical Neurons. Journal of Neurophysiology (2007) vol. 97 (3) pp. 2533 */
  double result = 0.0;

	// We are looking at the information flow of array1 ("I") -> array2 ("J")
	
	// clear memory
	memset(F_Ipast, 0, DATA_BINS*sizeof(unsigned long));
	for (char x=0; x<DATA_BINS; x++)
	{
		memset(F_Inow_Ipast[x], 0, DATA_BINS*sizeof(unsigned long));
		memset(F_Ipast_Jpast[x], 0, DATA_BINS*sizeof(unsigned long));
		for (char x2=0; x2<DATA_BINS; x2++)
			memset(F_Inow_Ipast_Jpast[x][x2], 0, DATA_BINS*sizeof(unsigned long));
	}
	
  // extract probabilities (actually number of occurrence)
	for (unsigned long t=WORD_LENGTH; t<NUM_SAMPLES; t++)
  {
		F_Ipast[arrayI[t-1]]++;
		F_Inow_Ipast[arrayI[t]][arrayI[t-1]]++;
		F_Ipast_Jpast[arrayI[t-1]][arrayJ[t-1]]++;
		F_Inow_Ipast_Jpast[arrayI[t]][arrayI[t-1]][arrayJ[t-1]]++;
#ifdef SHOW_DETAILED_PROGRESS
		status(t, REPORTS, NUM_SAMPLES-WORD_LENGTH);
#endif
	}
	
	for (char k=0; k<DATA_BINS; k++)
		for (char l=0; l<DATA_BINS; l++)
		{
			if (F_Ipast_Jpast[k][l] > 0)
				for (char m=0; m<DATA_BINS; m++)
					// if (F_Ipast[m]*F_Inow_Ipast[m][l]*F_Inow_Ipast_Jpast[m][k][l] != 0)
					if (F_Inow_Ipast_Jpast[m][k][l] != 0)
						result += F_Inow_Ipast_Jpast[m][k][l]/double(NUM_SAMPLES-WORD_LENGTH) * \
							log(double(F_Inow_Ipast_Jpast[m][k][l]*F_Ipast[k])/(F_Ipast_Jpast[k][l]*F_Inow_Ipast[m][k]));
		}

  return result/log(2);
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
