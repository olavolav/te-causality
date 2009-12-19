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

#define REPORTS 25
#define SHOW_DETAILED_PROGRESS
#define INFORMATION_IN_BITS

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
#define NUM_NEURONS 100
#define NUM_SAMPLES 126732
#define DATA_BINS 15
#define INPUTFILE "output/xresponse_10ms_15bins.dat"
#define OUTPUTFILE "output/transferentropy_os_10ms_15bins.mx"
// DemianTest with 20ms sampling
// #define NUM_NEURONS 100
// #define NUM_SAMPLES 89800
// #define DATA_BINS 5
// #define INPUTFILE "output/xresponse_5bins.dat"
// #define OUTPUTFILE "te/transferentropy_os_5bins.mx"

using namespace std;

double TEterm(char *array1, char *array2, char k, char l, char m);
double TransferEntropy(char *array1, char *array2);
void write_result(double **array);
void load_data(char **array);

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
    memset(xdata[i], 0, (NUM_SAMPLES)*sizeof(char));
    xresult[i] = new double[NUM_NEURONS];
    memset(xresult[i], 0, (NUM_NEURONS)*sizeof(double));
  }
  cout <<" done."<<endl;

  cout <<"loading data..."<<flush;
  load_data(xdata);
  cout <<" done."<<endl;

  // main loop:
  totaltrials = NUM_NEURONS*(NUM_NEURONS-1);
  cout  <<"set-up: "<<NUM_NEURONS<<" neurons, "<<totaltrials<<" trials"<<endl;
  completedtrials = 0;
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

  write_result(xresult);

  for(int i=0; i<NUM_NEURONS; i++) delete[] xdata[i];
  delete[] xdata;
  for(int i=0; i<NUM_NEURONS; i++) delete[] xresult[i];
  delete[] xresult;

  return 0;
}

double TEterm(char *array1, char *array2, char k, char l, char m)
{
  unsigned long countA, countB, countC, countD;
  countA = countB = countC = countD = 0;
  double result = 0.0;

  for (unsigned long tt=1; tt<NUM_SAMPLES; tt++)
  {
  	if (array2[tt-1] == l)
  	{
  		countD++;
  		if(array2[tt] == k) countC++;
  		if(array1[tt-1] == m)
  		{
  			countB++;
  			if(array2[tt] == k) countA++;
  		}
  	}
  }

  if (countA*countB*countC*countD != 0)
  {
  	result = (double)(countA)/NUM_SAMPLES * log((double)(countA*countD)/(countB*countC));
#ifdef INFORMATION_IN_BITS
  	// transform to information in bits
  	result /= log(2);
#endif
  }
  // else cout <<"!";

	return result;
}

double TransferEntropy(char *array1, char *array2)
{
	/* see for reference:
	     Gourevitch und Eggermont. Evaluating Information Transfer Between Auditory
	     Cortical Neurons. Journal of Neurophysiology (2007) vol. 97 (3) pp. 2533 */
  double result = 0.0;

  for (char k=0; k<DATA_BINS; k++)
  	for (char l=0; l<DATA_BINS; l++)
    {
#ifdef SHOW_DETAILED_PROGRESS
  		status(k*DATA_BINS+l, REPORTS, DATA_BINS*DATA_BINS);
#endif
      for (char m=0; m<DATA_BINS; m++)
      	result += TEterm(array1, array2, k, l, m);
    }

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
