// Created by olav, Di 30 Jun 2009 23:01:22 JST

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <ctime>
#include <cstring>
#include <sstream>

#include "../../Simulationen/olav.h"

#define REPORTS 25
#define SHOW_PROGRESS_BARS

#define NUM_NEURONS 100
#define NUM_SAMPLES 59951
#define DATA_BINS 30

using namespace std;

double TEterm(char *array1, char *array2, char k, char l, char m);
double TransferEntropy(char *array1, char *array2);
void write_result(double **array, string path);
void load_data(char **array);
void load_global_data(char *array);

int main(int argc, char *argv[])
{
  cout <<"------ transferentropy:te-global ------ olav, Di 30 Jun 2009 ------"<<endl;
  time_t start, end, middle;
  unsigned long totaltrials, completedtrials;

  time(&start);
  cout <<"start: "<<ctime(&start)<<flush;
  cout <<"running on host: "<<flush;
  system("hostname");
  cout <<"current directory: "<<flush;
  system("pwd");
  string opath="output/";
  cout <<"output directory: "<<opath<<endl;

  cout <<"allocating memory..."<<flush;
  char *xglobal = new char[NUM_SAMPLES];
  memset(xglobal, 0, (NUM_SAMPLES)*sizeof(char));
  char **xdata = new char*[NUM_NEURONS];
  double **xresult = new double*[NUM_NEURONS];
  for(int i=0; i<NUM_NEURONS; i++)
  {
    xdata[i] = new char[NUM_SAMPLES];
    memset(xdata[i], 0, (NUM_SAMPLES)*sizeof(char));
    xresult[i] = new double[2];
    memset(xresult[i], 0, 2*sizeof(double));
  }
  cout <<" done."<<endl;

  cout <<"loading data..."<<flush;
  load_data(xdata);
  load_global_data(xglobal);
  cout <<" done."<<endl;

  // main loop:
  totaltrials = 2*NUM_NEURONS;
  cout  <<"set-up: "<<NUM_NEURONS<<" neurons, "<<totaltrials<<" trials"<<endl;
  completedtrials = 0;

  for(int ii=0; ii<NUM_NEURONS; ii++)
	{
		cout <<"#"<<ii+1<<" -> global: "<<flush;
		xresult[ii][0] = TransferEntropy(xdata[ii], xglobal);
		time(&middle);
		cout <<" (elapsed "<<sec2string(difftime(middle,start))<<",";
		cout <<" ETA "<<sec2string((totaltrials-(completedtrials+1))*difftime(middle,start)/(completedtrials+1))<<")"<<endl;
		completedtrials++;

		cout <<"global -> #"<<ii+1<<": "<<flush;
		xresult[ii][1] = TransferEntropy(xglobal, xdata[ii]);
		time(&middle);
		cout <<" (elapsed "<<sec2string(difftime(middle,start))<<",";
		cout <<" ETA "<<sec2string((totaltrials-(completedtrials+1))*difftime(middle,start)/(completedtrials+1))<<")"<<endl;
		completedtrials++;
	}

  time(&end);
  cout <<"end: "<<ctime(&end)<<flush;
  cout <<"runtime: "<<sec2string(difftime(end,start))<<endl;

  write_result(xresult, opath);

  for(int i=0; i<NUM_NEURONS; i++) delete[] xdata[i];
  delete[] xdata;
  for(int i=0; i<2; i++) delete[] xresult[i];
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
  	result = (double)(countA)/NUM_SAMPLES * log((double)(countA*countD)/(countB*countC));

	return result;
}

double TransferEntropy(char *array1, char *array2)
{
  double result = 0.0;

  for (char k=0; k<DATA_BINS; k++)
  	for (char l=0; l<DATA_BINS; l++)
    {
      status(k*DATA_BINS+l, REPORTS, DATA_BINS*DATA_BINS);
      for (char m=0; m<DATA_BINS; m++)
      	result += TEterm(array1, array2, k, l, m);
    }

  return result;
}

void write_result(double **array, string path)
{
  ofstream fileout1((path+"transferentropy_global.mx").c_str());

  bool already;
  fileout1 <<"{";
  for(int j=0; j<NUM_NEURONS; j++)
  {
  	if(j>0) fileout1<<",";
  	fileout1 <<"{";
  	already = false;
    for(unsigned long i=0; i<2; i++)
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
  ifstream binaryfile("input/xresponseD.dat", ios::binary);
  if (binaryfile == NULL)
  {
  	cout <<"error: cannot find LOCAL input file!"<<endl;
  	exit(1);
  }

  for(int j=0; j<NUM_NEURONS; j++)
  {
    // status(j, REPORTS, NUM_NEURONS);
    binaryfile.read(array[j], NUM_SAMPLES);
  }
}

void load_global_data(char *array)
{
  ifstream binaryfile("input/xresponseglobalD.dat", ios::binary);
  if (binaryfile == NULL)
  {
  	cout <<"error: cannot find GLOBAL input file!"<<endl;
  	exit(1);
  }

  binaryfile.read(array, NUM_SAMPLES);
}
