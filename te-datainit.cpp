// created Do 28 Apr 2011 18:13:44 CEST by olav

#include "te-datainit.h"

#define FMODEL_SPIKECOUNT 1
#define FMODEL_HOWMANYAREACTIVE 2
#define FMODEL_LEOGANG 3
#define FMODEL_ERROR -1

#define HEIGHT_OF_ASCII_HISTOGRAMS 12


// set output stream depending on wether SimKernel's sim.h is included
// (see also te-datainit.h)
#undef IOSTREAMH
#undef IOSTREAMC
#undef IOSTREAMV

#ifdef SIM_IO_H
  // SimKernel found.
  #define IOSTREAMH Sim& output
  #define IOSTREAMC output.io
  #define IOSTREAMV output
  #define IOSTREAMENDL Endl
#else
  // SimKernel not found, using standard output.
  #define IOSTREAMH std::ostream* output
  #define IOSTREAMC *output
  #define IOSTREAMV output
  #define IOSTREAMENDL std::endl
#endif

// TEST
// #define IOSTREAMH std::ostream* output
// #define IOSTREAMC *output
// #define IOSTREAMV output


double** load_time_series_from_binary_file(std::string inputfile_name, unsigned int size, long samples, double input_scaling, bool OverrideRescalingQ, double std_noise, double fluorescence_saturation, double cutoff, IOSTREAMH)
{
	// initialize random number generator
	gsl_rng* GSLrandom;
	gsl_rng_env_setup();
	GSLrandom = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(GSLrandom, 1234);
	
	// reserve and clear memory for result ("try&catch" is still missing!)
  double **xresult = NULL;
  xresult = new double*[size];
  for(unsigned int i=0; i<size; i++)
  {
    xresult[i] = NULL;
    xresult[i] = new double[samples];
    memset(xresult[i], 0, samples*sizeof(double));
  }
	// assert((BytesPerDataPoint==1)||(BytesPerDataPoint==2)); //so far
  char* in_from_file_array = new char[samples];
	double* tempdoublearray = new double[samples];
  // memset(tempdoublearray, 0, samples*sizeof(double));
	// double xtemp;
  
	// open input file
	char* name = new char[inputfile_name.length()+1];
	strcpy(name,inputfile_name.c_str());
  std::ifstream binaryfile(name, std::ios::binary);
	delete[] name;

  if (binaryfile == NULL)
  {
  	IOSTREAMC <<IOSTREAMENDL<<"error: cannot find input file!"<<IOSTREAMENDL;
  	exit(1);
  }

	// test file length
	binaryfile.seekg(0,std::ios::end);
	if(long(binaryfile.tellg()) != size*samples)
	{
  	IOSTREAMC <<IOSTREAMENDL<<"error: file length of input does not match given parameters!"<<IOSTREAMENDL;
  	exit(1);
	}
	binaryfile.seekg(0,std::ios::beg);
	
  for(int j=0; j<size; j++)
  {
    binaryfile.read(in_from_file_array, samples);

		// OverrideRescalingQ = true
		// Dies ignoriert also "appliedscaling", "noise", "HighPassFilterQ" und "cutoff"
		// Therefore, "bins" takes the role of an upper cutoff
		if (OverrideRescalingQ)
	    for(long k=0; k<samples; k++)
        // xdata[j][k] = in_from_file_array[k];
			  xresult[j][k] = double(in_from_file_array[k]);
		else     // OverrideRescalingQ = false
		{
	    for (long k=0; k<samples; k++)
			{
				// transform to unsigned notation
				tempdoublearray[k] = double(in_from_file_array[k]);
				if (in_from_file_array[k]<0) tempdoublearray[k] += 256.;
				
				// transform back to original signal (same as in Granger case)
				tempdoublearray[k] /= input_scaling;
				// assuming a saturation with hill function of order 1
				if (fluorescence_saturation > 0.)
					tempdoublearray[k] = tempdoublearray[k]/(tempdoublearray[k] + fluorescence_saturation);
				// adding noise
				if (std_noise > 0.)
					tempdoublearray[k] += gsl_ran_gaussian(GSLrandom,std_noise);					
				// apply cutoff
				if ((cutoff>0)&&(tempdoublearray[k]>cutoff)) tempdoublearray[k] = cutoff;
			}
			
		}
    memcpy(xresult[j],tempdoublearray,samples*sizeof(double));
  }


  // // determine available samples per globalbin for TE normalization later
  // memset(AvailableSamples, 0, globalbins*sizeof(unsigned long));
  // for (unsigned long t=StartSampleIndex; t<EndSampleIndex; t++)
  //  AvailableSamples[xglobal[t]]++;
  // 
  // if (EqualSampleNumberQ || (MaxSampleNumberPerBin>0))
  // {
  //  unsigned long maxsamples = ULONG_MAX;
  //  for (rawdata g=0; g<globalbins; g++)
  //    if (AvailableSamples[g]<maxsamples) maxsamples = AvailableSamples[g];
  //  IOSTREAMC <<"DEBUG: maxsamples = "<<maxsamples<<IOSTREAMENDL;
  //  
  //  if ((MaxSampleNumberPerBin>maxsamples) && !EqualSampleNumberQ)
  //    maxsamples = MaxSampleNumberPerBin;
  //  if ((MaxSampleNumberPerBin<maxsamples)&&(MaxSampleNumberPerBin>0))
  //    maxsamples = MaxSampleNumberPerBin;
  //  IOSTREAMC <<"DEBUG: cut to maxsamples = "<<maxsamples<<IOSTREAMENDL;
  //  
  //  unsigned long* AlreadySelectedSamples = new unsigned long[globalbins];
  //  memset(AlreadySelectedSamples, 0, globalbins*sizeof(unsigned long));
  //  for (unsigned long t=StartSampleIndex; t<EndSampleIndex; t++)
  //    if ((++AlreadySelectedSamples[xglobal[t]])>maxsamples)
  //      xglobal[t] = globalbins; // ..and therefore exclude from calculation
  // 
  //  // re-determine available samples per globalbin (inefficient)
  //  memset(AvailableSamples, 0, globalbins*sizeof(unsigned long));
  //  for (unsigned long t=StartSampleIndex; t<EndSampleIndex; t++)
  //    if (xglobal[t]<globalbins) AvailableSamples[xglobal[t]]++;
  //    
  //  delete[] AlreadySelectedSamples;
  // }
	
	// free allocated memory
	gsl_rng_free(GSLrandom);    
  delete[] in_from_file_array;
  // delete[] xglobaltemp;
	delete[] tempdoublearray;
  // delete[] tempdoublearraycopy;
	
  return xresult;
};


rawdata* generate_discretized_global_time_series(double** time_series, unsigned int size, long samples, unsigned int globalbins, double GlobalConditioningLevel, unsigned long* AvailableSamples, long StartSampleIndex, long EndSampleIndex, IOSTREAMH)
{
	rawdata* xglobal = new rawdata[samples];
	memset(xglobal, 0, samples*sizeof(rawdata));
  double* xglobaltemp = generate_mean_time_series(time_series, size, samples);
  	
	// EVIL SAALBACH HACK FOR TIME CODE GLOBAL SIGNAL: -------------------------------------------- !!!!!!!!!!
	// xglobaltemp[0] = 0.;
	// for (unsigned long t=0; t<samples; t++)
	// 	xglobaltemp[t] = double(int(t)%int(60*24/tauF));
	// 	// xglobaltemp[t] = double(mod(t,60*24/tauF));
	// IOSTREAMC <<"DEBUG OF EVIL TIME CODE HACK: last globaltemp value = "<<xglobaltemp[samples-1]<<IOSTREAMENDL;
	
	if (GlobalConditioningLevel > 0.)
	{
		unsigned long below = 0;
		for (long t=0; t<samples; t++)
		{
			if (xglobaltemp[t] > GlobalConditioningLevel) xglobal[t] = 1;
			else
			{
				xglobal[t] = 0;
				below++;
			}
		}
    IOSTREAMC <<" -> global conditioning level "<<GlobalConditioningLevel<<": "<<(100.*below)/samples;
    IOSTREAMC <<"% are below threshold. "<<IOSTREAMENDL;
	}
	else discretize(xglobaltemp,xglobal,samples,globalbins);
	
	// determine available samples per globalbin for TE normalization later
  memset(AvailableSamples, 0, globalbins*sizeof(unsigned long));
  for (unsigned long t=StartSampleIndex; t<EndSampleIndex; t++)
   AvailableSamples[xglobal[t]]++;


  // if (EqualSampleNumberQ || (MaxSampleNumberPerBin>0))
  // {
  //  unsigned long maxsamples = ULONG_MAX;
  //  for (rawdata g=0; g<globalbins; g++)
  //    if (AvailableSamples[g]<maxsamples) maxsamples = AvailableSamples[g];
  //  IOSTREAMC <<"DEBUG: maxsamples = "<<maxsamples<<IOSTREAMENDL;
  //  
  //  if ((MaxSampleNumberPerBin>maxsamples) && !EqualSampleNumberQ)
  //    maxsamples = MaxSampleNumberPerBin;
  //  if ((MaxSampleNumberPerBin<maxsamples)&&(MaxSampleNumberPerBin>0))
  //    maxsamples = MaxSampleNumberPerBin;
  //  IOSTREAMC <<"DEBUG: cut to maxsamples = "<<maxsamples<<IOSTREAMENDL;
  //  
  //  unsigned long* AlreadySelectedSamples = new unsigned long[globalbins];
  //  memset(AlreadySelectedSamples, 0, globalbins*sizeof(unsigned long));
  //  for (unsigned long t=StartSampleIndex; t<EndSampleIndex; t++)
  //    if ((++AlreadySelectedSamples[xglobal[t]])>maxsamples)
  //      xglobal[t] = globalbins; // ..and therefore exclude from calculation
  // 
  //  // re-determine available samples per globalbin (inefficient)
  //  memset(AvailableSamples, 0, globalbins*sizeof(unsigned long));
  //  for (unsigned long t=StartSampleIndex; t<EndSampleIndex; t++)
  //    if (xglobal[t]<globalbins) AvailableSamples[xglobal[t]]++;
  //    
  //  delete[] AlreadySelectedSamples;
  // }


  free_time_series_memory(xglobaltemp);
  return xglobal;
}

rawdata** generate_discretized_version_of_time_series(double** in, unsigned int size, long nr_samples, unsigned int nr_bins)
{
  rawdata** xout;
  xout = new rawdata*[size];
  for(unsigned int ii=0; ii<size; ii++)
  {
    xout[ii] = new rawdata[nr_samples];
    discretize(in[ii], xout[ii], nr_samples, nr_bins);
  }  
  return xout;
};

void discretize(double* in, rawdata* out, long nr_samples, unsigned int nr_bins)
{
  discretize(in,out,smallest(in,nr_samples),largest(in,nr_samples),nr_samples,nr_bins);
};
void discretize(double* in, rawdata* out, double min, double max, long nr_samples, unsigned int nr_bins)
{
	double xstepsize = (max-min)/nr_bins;

	int xint;
	for (unsigned long t=0; t<nr_samples; t++)
    out[t] = discretize(in[t],min,max,nr_bins);
};
rawdata discretize(double in, double min, double max, unsigned int nr_bins)
{
  assert(max>min);
  assert(nr_bins>0);
  // correct discretization according to 'te-test.nb'
  // incorporated later: double xstepsize = (max-min)/nr_bins;
  // IOSTREAMC <<"max = "<<max<<IOSTREAMENDL;
  // IOSTREAMC <<"min = "<<min<<IOSTREAMENDL;
  // IOSTREAMC <<"bins here = "<<nr_bins<<IOSTREAMENDL;
  // IOSTREAMC <<"stepsize = "<<xstepsize<<IOSTREAMENDL;

  int xint;

  // assert(in<=max); ...does not have to be true, and does not matter, data is included in highest bin then
  // assert(in>=min);
	if (in>=max) xint = nr_bins-1;
	else
	{
		if (in<=min) xint = 0;
		// with stepsize variable: else xint = (int)((in-min)/xstepsize);
		// without:
		else xint = (int)((in-min)*double(nr_bins)/(max-min));
	}
	if (xint >= nr_bins) xint = nr_bins-1; // need to have this for some silly numerical reason...

	assert((xint>=0)&&(rawdata(xint)<nr_bins)); // ...just to be sure...
	return rawdata(xint);
};

// #ifdef ENABLE_ADAPTIVE_BINNING_AT_COMPILE_TIME 
// void discretize2accordingtoStd(double* in, rawdata* out)
// {
//  double splitheight = sqrt(gsl_stats_variance(in,1,samples));
// 
//  int xint;
//  for (unsigned long t=0; t<samples; t++)
//  {
//    if (in[t]>splitheight) out[t] = 1;
//    else out[t] = 0;
//  }
// };
// #endif

void apply_high_pass_filter_to_time_series(double** time_series, unsigned int size, long nr_samples)
{
  for(unsigned int ii=0; ii<size; ii++)
    apply_high_pass_filter_to_time_series(time_series[ii], nr_samples);
};
void apply_high_pass_filter_to_time_series(double* time_series, long nr_samples)
{
	double* arraycopy = new double[nr_samples];

  // of course, this is just a difference signal, so not really filtered
  memcpy(arraycopy,time_series,nr_samples*sizeof(double));
  time_series[0] = 0.;
  for(long k=1; k<nr_samples; k++)
  	time_series[k] = arraycopy[k] - arraycopy[k-1];

  delete[] arraycopy;
};

double** generate_time_series_from_spike_data(std::string inputfile_spiketimes, std::string inputfile_spikeindices, unsigned int size, unsigned int tauImg, long samples, std::string fluorescence_model, double std_noise, double fluorescence_saturation, double cutoff, double DeltaCalciumOnAP, double tauCa)
{
	// reserve and clear memory for result ("try&catch" is still missing!)
  double **xresult = new double*[size];
  for(unsigned int i=0; i<size; i++)
  {
    xresult[i] = new double[samples];
    memset(xresult[i], 0, samples*sizeof(double));
  }
  
	// open files
	char* nameI = new char[inputfile_spikeindices.length()+1];
	strcpy(nameI,inputfile_spikeindices.c_str());
  std::ifstream binaryfileI(nameI, std::ios::binary);
	delete[] nameI;
	char* nameT = new char[inputfile_spiketimes.length()+1];
	strcpy(nameT,inputfile_spiketimes.c_str());
  std::ifstream binaryfileT(nameT, std::ios::binary);
	delete[] nameT;
	
	// determine file length
	binaryfileI.seekg(0,std::ios::end);
	const long nr_spikes = binaryfileI.tellg()/sizeof(int);
  // IOSTREAMC <<"number of spikes in index file: "<<nr_spikes<<IOSTREAMENDL;
	binaryfileI.seekg(0,std::ios::beg);
	int* xindex = new int[nr_spikes];
	double* xtimes = new double[nr_spikes];
	
	// read spike data
	binaryfileI.read((char*)xindex, nr_spikes*sizeof(int));
	binaryfileT.read(reinterpret_cast<char*>(xtimes), nr_spikes*sizeof(double));
	
	// close files
  binaryfileI.close();
  binaryfileT.close();
	
	// test if read data appears valid
  for (long tt=0; tt<samples; tt++)
  {
    assert((xindex[tt]>=0)&&(xindex[tt]<size)); // indices are in allowed range
    if(tt>0) assert(xtimes[tt]>=xtimes[tt-1]); // spike times are an ordered sequence
    if(tt<samples-1) assert(xtimes[tt]<=xtimes[tt+1]);
  }
	
  // for(long t=0; t<20; t++)
  //   IOSTREAMC <<"DEBUG: xindex = "<<xindex[t]<<", xtimes = "<<xtimes[t]<<IOSTREAMENDL;

  // choose switch key for the fluorescence model
  int fluorescence_model_key = FMODEL_ERROR;
  if (fluorescence_model == "SpikeCount") fluorescence_model_key = FMODEL_SPIKECOUNT;
  if (fluorescence_model == "HowManyAreActive") fluorescence_model_key = FMODEL_HOWMANYAREACTIVE;
  if (fluorescence_model == "Leogang") fluorescence_model_key = FMODEL_LEOGANG;
  assert(fluorescence_model_key != FMODEL_ERROR);
  
	// generate fluorescence data
  // const int int_tauF = (int)round(tauF); // in ms
  unsigned long startindex = 1;
  unsigned long endindex = 0; // therefore, we miss the first spike!
  long dataindex = 0;
  unsigned long tinybit_spikenumber;
  // unsigned long ttExactMS = 0;
  for (unsigned long ttExactMS=0; ttExactMS<tauImg*samples; ttExactMS+=tauImg)
  {
    // determine starting and ending spike index of current frame
    while ((endindex+1<nr_spikes)&&(xtimes[endindex+1]<=ttExactMS+tauImg))
      endindex++;
    tinybit_spikenumber = std::max((int)(endindex-startindex+1),0);
    
    // IOSTREAMC <<"DEBUG: ttExactMS = "<<ttExactMS<<", startindex = "<<startindex<< \
    //   ", endindex = "<<endindex<<", tinybit_spikenumber = "<<tinybit_spikenumber<<IOSTREAMENDL;

    for (int ii=0; ii<size; ii++)
    {
      switch (fluorescence_model_key)
      {
        case FMODEL_SPIKECOUNT:
          if(tinybit_spikenumber>0)
            xresult[ii][dataindex] = double(count(xindex,startindex,endindex,ii));
            // test: xresult[ii][dataindex] = double(tinybit_spikenumber);
          break;
          
        case FMODEL_HOWMANYAREACTIVE:
          if(tinybit_spikenumber>0)
            xresult[ii][dataindex] = double(has_index(xindex,startindex,endindex,ii));
          break;
          
        case FMODEL_LEOGANG:
          xresult[ii][dataindex] = (1.-double(tauImg)/tauCa)*xresult[ii][std::max(dataindex-1,long(0))] + \
            DeltaCalciumOnAP*double(count(xindex,startindex,endindex,ii));
          break;
          
        // default:
        //   IOSTREAMC <<"error in generate_time_series_from_spike_data: invalid fluorescence model"<<IOSTREAMENDL;
        //   exit(1);
      }
    }
    if(startindex <= endindex)
      startindex = 1 + endindex;
    dataindex++;
  }

  // apply saturation (Hill function of order 1 as usual)
  if(fluorescence_saturation > 0.)
    for (unsigned int ii=0; ii<size; ii++)
      for (long tt=0; tt<samples; tt++)
        xresult[ii][tt] = xresult[ii][tt]/(xresult[ii][tt]+fluorescence_saturation);

  // apply additive noise term
  if(std_noise > 0.)
  {
    // initialize random number generator
  	gsl_rng* GSLrandom;
  	gsl_rng_env_setup();
  	GSLrandom = gsl_rng_alloc(gsl_rng_default);
  	gsl_rng_set(GSLrandom, 1234);
  	
    for (unsigned int ii=0; ii<size; ii++)
      for (long tt=0; tt<samples; tt++)
        xresult[ii][tt] += gsl_ran_gaussian(GSLrandom,std_noise);
        
		// free allocated memory
		gsl_rng_free(GSLrandom);    
  }
  
	delete[] xindex;
	delete[] xtimes;
	
  return xresult;
};

unsigned long count(int* array, unsigned long starti, unsigned long endi, int what)
{
	unsigned long occur = 0;
	for (unsigned long i=starti; i<=endi; i++)
		if (array[i] == what) occur++;
	return occur;
};
bool has_index(int* array, unsigned long starti, unsigned long endi, int what)
{
	for (unsigned long i=starti; i<=endi; i++)
		if (array[i] == what) return true;
	return false;
};

double smallest(double* array, const long length)
{
	double min = array[0];
	for (long i=1; i<length; i++)
		if(array[i]<min) min = array[i];

	return min;
};
double largest(double* array, const long length)
{
	double max = array[0];
	for (long i=1; i<length; i++)
		if(array[i]>max) max = array[i];

	return max;
};
rawdata smallest(rawdata* array, const long length)
{
	rawdata min = array[0];
	for (long i=1; i<length; i++)
		if(array[i]<min) min = array[i];

	return min;
};
rawdata largest(rawdata* array, const long length)
{
	rawdata max = array[0];
	for (long i=1; i<length; i++)
		if(array[i]>max) max = array[i];

	return max;
};

double smallest(double** array, const unsigned int size, const long length)
{
	double min = array[0][0];
	for (unsigned int i=1; i<size; i++)
    min = std::min(min,smallest(array[i],length));

	return min;
};
double largest(double** array, const unsigned int size, const long length)
{
	double max = array[0][0];
	for (unsigned int i=1; i<size; i++)
    max = std::max(max,largest(array[i],length));

	return max;
};

double total(double* array, const long length)
{
	double sum = 0.;
	for (long i=1; i<length; i++)
		sum += array[i];

	return sum;
};

double* generate_mean_time_series(double** data, unsigned int size, long samples)
{
	double* xglobaltemp = new double[samples];
	memset(xglobaltemp, 0, samples*sizeof(double));

	for (long t=0; t<samples; t++)
  {
    for (unsigned int ii=0; ii<size; ii++)
  		xglobaltemp[t] += data[ii][t];
		xglobaltemp[t] /= size;
  }
  
  return xglobaltemp;
}

void free_time_series_memory(double** xresult, unsigned int size)
{
  for(unsigned int ii=0; ii<size; ii++)
    free_time_series_memory(xresult[ii]);
  delete[] xresult;
};
void free_time_series_memory(double* xresult)
{
  delete[] xresult;
};
void free_time_series_memory(rawdata** xresult, unsigned int size)
{
  for(unsigned int ii=0; ii<size; ii++)
    free_time_series_memory(xresult[ii]);
  delete[] xresult;
};
void free_time_series_memory(rawdata* xresult)
{
  delete[] xresult;
};

void display_subset(double* data, int length, IOSTREAMH)
{
  // IOSTREAMC <<"displaying some subset of data points:"<<IOSTREAMENDL;
  IOSTREAMC <<"{";
  for (long t=0; t<length; t++)
  {
    if (t>0) IOSTREAMC <<",";
    IOSTREAMC <<data[t];
  }
  IOSTREAMC <<"} (range "<<smallest(data,length)<<" – "<<largest(data,length)<<")"<<IOSTREAMENDL;
};
void display_subset(rawdata* data, int length, IOSTREAMH)
{
  // IOSTREAMC <<"displaying some subset of data points:"<<IOSTREAMENDL;
  IOSTREAMC <<"{";
  for (long t=0; t<length; t++)
  {
    if (t>0) IOSTREAMC <<",";
    IOSTREAMC <<int(data[t]);
  }
  IOSTREAMC <<"} (range "<<int(smallest(data,length))<<" – "<<int(largest(data,length))<<")"<<IOSTREAMENDL;
};

int Magic_GuessBinNumber(double** data, const unsigned int size, const long samples)
{
  double range, std;  
  double meanbins = 0.;
  for(unsigned int i = 0; i < size; i++)
  {
    range = largest(data[i],samples)-smallest(data[i],samples);
    assert(range > 0.);
    std += sqrt(gsl_stats_variance(data[i],1,samples));
    meanbins += 2*std/range;
  }
  meanbins /= size;
  
  return std::max(2,int(round(meanbins)));
};

double Magic_GuessConditioningLevel(double** data, const unsigned int size, const long samples, IOSTREAMH)
{
  int histo_bins = int(std::max(4.,std::min(200.,sqrt(samples))));
  IOSTREAMC <<" -> number of bins for histogram: "<<histo_bins<<IOSTREAMENDL;
  double xmeanmin, xmeanmax;
  double xresultlevel = -1.;
  
  double* xmean = generate_mean_time_series(data,size,samples);
  xmeanmin = smallest(xmean,samples);
  xmeanmax = largest(xmean,samples);
  // IOSTREAMC <<"-> xmeanmin = "<<xmeanmin<<", xmeanmax = "<<xmeanmax<<IOSTREAMENDL;
  // IOSTREAMC <<" -> beginning of <f>: "<<IOSTREAMENDL;
  // display_subset(xmean,5,IOSTREAMV);
  
  // find maximum, which we assume comes from the noise peak
  double x_ymax = Util_FindPeakInHistogram(xmean,samples,xmeanmin,xmeanmax,histo_bins);
  IOSTREAMC <<" -> identified peak at <f> = "<<x_ymax<<IOSTREAMENDL;
  
  // for(int i=0; i<histo_bins; i++)
  //   IOSTREAMC <<"histo <f> = "<<xmeanmin+(double(i)+0.5)*(xmeanmax-xmeanmin)/double(histo_bins)<<": count = "<<histo[i]<<IOSTREAMENDL;
  
  // PlotHistogramInASCII(xmean,samples,xmeanmin,x_ymax+0.1,"<f>","#(<f>)",IOSTREAMV);
  IOSTREAMC <<IOSTREAMENDL<<" -> log histogram of complete range of <f>:"<<IOSTREAMENDL;
  PlotLogHistogramInASCII(xmean,samples,xmeanmin,xmeanmax,"<f>","log #(<f>)",IOSTREAMV);
  // IOSTREAMC <<IOSTREAMENDL<<" -> log histogram of the right tail of <f>:"<<IOSTREAMENDL;
  // PlotLogHistogramInASCII(xmean,samples,x_ymax,xmeanmax,"log <f>","log #(<f>)",IOSTREAMV);
  IOSTREAMC <<IOSTREAMENDL<<" -> log-log histogram of the right tail of <f>:"<<IOSTREAMENDL;
  PlotLogLogHistogramInASCII(xmean,samples,x_ymax,xmeanmax,"log <f>","log #(<f>)",IOSTREAMV);
  
  // re-sample right tail of histogram (new method)
  double *x = NULL;
  double *y = NULL;
  double *w = NULL;
  double c0,c1,cov00,cov01,cov11,chisq;
  long nr_observations = long(0.2*round(sqrt(samples)));
  Util_CreateFakeLogLogHistogram(&x,&y,&w,xmean,samples,x_ymax,xmeanmax,nr_observations);
  
  // make linear weighted fit using GSL
  gsl_fit_wlinear(x,1,w,1,y,1,nr_observations,&c0,&c1,&cov00,&cov01,&cov11,&chisq);
  IOSTREAMC <<" -> best fit: y(x) = "<<c0<<" + ("<<c1<<")*x; (chisq = "<<chisq<<")"<<IOSTREAMENDL;

  // IOSTREAMC <<"vectorxy = ";
  // Util_CoordinatedForMathematica(x,y,nr_observations,IOSTREAMV);
  // IOSTREAMC <<"vectorxw = ";
  // Util_CoordinatedForMathematica(x,w,nr_observations,IOSTREAMV);
  // IOSTREAMC <<"c0 = "<<c0<<";"<<IOSTREAMENDL;
  // IOSTREAMC <<"c1 = "<<c1<<";"<<IOSTREAMENDL<<IOSTREAMENDL;
  
  // tranforming everything back to linear space (where the line in log-log corresponds
  // to the power-law fit we wanted)
  IOSTREAMC <<" -> transforming back to linear space:"<<IOSTREAMENDL;
  for(int i=0; i<nr_observations; i++)
  {
    x[i] = exp(x[i]);
    y[i] = exp(y[i]);
    // w[i] = exp(w[i]); Unsinn, das ist ja normiert und so.

    // set up new weights in linear space
    w[i] = exp(-pow((x[i]-(x_ymax+xmeanmax)/2.)/((xmeanmax-x_ymax)/4.),2.));
  }
  // normalize weights (actually, this would not be necessary but is nicer maybe.)
  w[0] = 0.; // evil hack for the peak... :-(
  for(int i=0; i<histo_bins; i++)
    w[i] = w[i]/total(w,nr_observations);
    
  double coeffA = exp(c0);
  double coeffGamma = c1;
  IOSTREAMC <<" -> best fit in linear space: p(f) = "<<coeffA<<" * f^("<<coeffGamma<<")"<<IOSTREAMENDL;

  double xthresh;
  double* deviations = new double[nr_observations];
  for(int i=0; i<nr_observations; i++)
    deviations[i] = coeffA*pow(x[i],coeffGamma)-y[i];
  xthresh = 0.5*sqrt(gsl_stats_wvariance(w,1,deviations,1,nr_observations));
  IOSTREAMC <<" -> threshold for deviation: "<<xthresh<<IOSTREAMENDL;

  // IOSTREAMC <<" -> x = ";
  // display_subset(x,nr_observations,IOSTREAMV);
  // IOSTREAMC <<" -> w = ";
  // display_subset(w,nr_observations,IOSTREAMV);
  // IOSTREAMC <<" -> deviations = ";
  // display_subset(deviations,nr_observations,IOSTREAMV);

  IOSTREAMC <<"vectorxy = ";
  Util_CoordinatedForMathematica(x,y,nr_observations,IOSTREAMV);
  IOSTREAMC <<"vectorxw = ";
  Util_CoordinatedForMathematica(x,w,nr_observations,IOSTREAMV);
  IOSTREAMC <<"A = "<<coeffA<<";"<<IOSTREAMENDL;
  IOSTREAMC <<"gamma = "<<coeffGamma<<";"<<IOSTREAMENDL<<IOSTREAMENDL;

  // Von der Hälfte der Werte (auf lin. Skala) ab gehen wir nach unten (f=0), bis die
  // Abweichung vom Fit überschwellig wird:
  bool FoundConditioningLevel = false;
  int debug_count = 0;
  if(coeffGamma<0.)
  {
    for(int i=0; i<nr_observations; i++)
    {
      // find greatest x on left side for which the deviations are above theshold
      if(((x[i]-x_ymax)/(xmeanmax-x_ymax)<0.5) && (abs(deviations[i])>xthresh))
      {
        debug_count++;
        if(!FoundConditioningLevel)
        {
          FoundConditioningLevel = true;
          xresultlevel = x[i];
        }
        else if(x[i]>xresultlevel)
          xresultlevel = x[i];
      }
      else break;
    }
  }
  else IOSTREAMC <<" -> Warning: power law exponent is not negative as expected!"<<IOSTREAMENDL;
  IOSTREAMC <<" -> Debug: found "<<debug_count<<" points with superthreshold deviations."<<IOSTREAMENDL;
  
  if(!FoundConditioningLevel)
  {
    IOSTREAMC <<" -> Warning: Conditioning level could not be found, override at 1/3 between max and end!";
    IOSTREAMC <<IOSTREAMENDL;
    xresultlevel = 0.33333*(xmeanmax-x_ymax) + x_ymax;
  }
  
  // IOSTREAMC <<" -> number of observations for variance estimate: "<<nr_observations<<IOSTREAMENDL;
  // IOSTREAMC <<" -> std. deviation of difference to power-law fit (pseudo-weighted): ";
  // IOSTREAMC <<sqrt(gsl_stats_variance(deviations,1,nr_observations))<<IOSTREAMENDL;
  
  // free memory
  Util_FreeFakeHistogram(x,y,w);
  free_time_series_memory(xmean);
  delete[] deviations;
  return xresultlevel;
};

void Test_SetMinimaToZero(double** data, unsigned int size, long samples)
{
  double minimum;
  for(unsigned int i=0; i<size; i++)
  {
    minimum = smallest(data[i],samples);
    for(long t=0; t<samples; t++) data[i][t] -= minimum;
    assert(!(smallest(data[i],samples)<0.));
  }
};

// void Test_PlotLogHistogram(long* histo, int length)
// {
//   long* loghisto = new long[length];
//   for(int i=0; i<length; i++)
//   {
//     if(histo[i]>3) loghisto[i] = long(round(log(histo[i])));
//     else loghisto[i] = 0;
//   }
//   Test_PlotHistogram(loghisto,length);
//   delete[] loghisto;
// };
// void Test_PlotHistogram(long* histo, int length)
// {
//  long max = histo[0];
//  for (unsigned int i=1; i<length; i++)
//     max = std::max(max,histo[i]);
//     
//   IOSTREAMC <<"^"<<IOSTREAMENDL;
//   for(int line=0; line<HEIGHT_OF_ASCII_HISTOGRAMS; line++)
//   {
//     IOSTREAMC <<"|";
//     for(int row=0; row<length; row++)
//     {
//       if(double(histo[row])/double(max)>(1.-double(line)/double(HEIGHT_OF_ASCII_HISTOGRAMS))) IOSTREAMC <<"#";
//       else IOSTREAMC <<" ";
//     }
//     IOSTREAMC <<IOSTREAMENDL;
//   }
//   IOSTREAMC <<"+";
//   for(int row=0; row<length; row++) IOSTREAMC <<"-";
//   IOSTREAMC <<">"<<IOSTREAMENDL;
// };

void PlotHistogramInASCII(double* data, int samples, double range_min, double range_max, std::string xlabel, std::string ylabel, IOSTREAMH)
{
  PlotHistogramInASCII(false,false,data,samples,range_min,range_max,xlabel,ylabel,IOSTREAMV);
};
void PlotLogHistogramInASCII(double* data, int samples, double range_min, double range_max, std::string xlabel, std::string ylabel, IOSTREAMH)
{
  PlotHistogramInASCII(false,true,data,samples,range_min,range_max,xlabel,ylabel,IOSTREAMV);
};
void PlotLogLogHistogramInASCII(double* data, int samples, double range_min, double range_max, std::string xlabel, std::string ylabel, IOSTREAMH)
{
  PlotHistogramInASCII(true,true,data,samples,range_min,range_max,xlabel,ylabel,IOSTREAMV);
};
void PlotHistogramInASCII(bool xlogscaleQ, bool ylogscaleQ, double* data, int samples, double range_min, double range_max, std::string xlabel, std::string ylabel, IOSTREAMH)
{
  assert(range_min!=range_max);
  const int histo_bins = std::max(20,std::min(65,int(round(sqrt(samples)))));
  // IOSTREAMC <<" -> debug: histo_bins = "<<histo_bins<<IOSTREAMENDL;
  
  // create and fill histogram
  long* histo = new long[histo_bins];
  memset(histo, 0, histo_bins*sizeof(long));
  if(!xlogscaleQ)
    for(long t=0; t<samples; t++)
      histo[int(discretize(data[t],range_min,range_max,histo_bins))] += 1;
  else
    for(long t=0; t<samples; t++)
      histo[int(discretize(log(data[t]),log(range_min),log(range_max),histo_bins))] += 1;

  // convert histogram to floating point array (and free integer histogram)
  double* histoD = new double[histo_bins];
  for (int i=0; i<histo_bins; i++)
    histoD[i] = double(histo[i]);
  delete[] histo;

  if(ylogscaleQ)
    for(int i=0; i<histo_bins; i++)
    {
      if(histoD[i]>exp(1.)) histoD[i] = log(histoD[i]);
      else histoD[i] = 0.;
    }
  
  // find maximum of histogram
  double max_histo = 0.;
	for (unsigned int i=0; i<histo_bins; i++)
    if(histoD[i]>max_histo) max_histo = histoD[i];
  long min_histo = max_histo;
	for (unsigned int i=0; i<histo_bins; i++)
    if(histoD[i]<min_histo) min_histo = histoD[i];
       
  // draw histogram
  IOSTREAMC <<"^";
  if(ylogscaleQ) IOSTREAMC <<" (log)";
  IOSTREAMC <<IOSTREAMENDL;
  for(int line=0; line<HEIGHT_OF_ASCII_HISTOGRAMS; line++)
  {
    IOSTREAMC <<"|";
    for(int row=0; row<histo_bins; row++)
    {
      // if(histoD[row]/double(max_histo)>=(1.-double(line)/double(HEIGHT_OF_ASCII_HISTOGRAMS))) IOSTREAMC <<"#";
      if(histoD[row]/max_histo>=(1.-(double(line)+0.5)/double(HEIGHT_OF_ASCII_HISTOGRAMS)))
      {
        if(histoD[row]/max_histo>=(1.-(double(line))/double(HEIGHT_OF_ASCII_HISTOGRAMS)))
          IOSTREAMC <<":";
        else IOSTREAMC <<".";
      }
      else IOSTREAMC <<" ";
    }
    IOSTREAMC <<IOSTREAMENDL;
  }
  IOSTREAMC <<"+";
  for(int row=0; row<histo_bins; row++) IOSTREAMC <<"-";
  IOSTREAMC <<">";
  if(xlogscaleQ) IOSTREAMC <<" (log)";
  IOSTREAMC <<IOSTREAMENDL;
  
  // label axes
  IOSTREAMC <<"x-axis: "<<xlabel<<", ";
  IOSTREAMC <<"range "<<range_min<<" – "<<range_max<<IOSTREAMENDL;
  IOSTREAMC <<"y-axis: "<<ylabel<<", ";
  if(!ylogscaleQ) IOSTREAMC <<"range "<<min_histo<<" – "<<max_histo<<IOSTREAMENDL;
  else IOSTREAMC <<"range "<<long(exp(min_histo))<<" – "<<long(exp(max_histo))<<IOSTREAMENDL;
  
  // free memory
  delete[] histoD;
};

double Util_FindPeakInHistogram(const double* data, const long samples, const double range_min, const double range_max, const int histo_bins)
{
  // create and fill histogram
  long* histo = new long[histo_bins];
  memset(histo, 0, histo_bins*sizeof(long));
  for(long t=0; t<samples; t++)
    histo[int(discretize(data[t],range_min,range_max,histo_bins))] += 1;

  double x = range_min;
  long y = 0;
  
  for(int i=0; i<histo_bins; i++)
  {
    if(histo[i]>y)
    {
      // found a higher point:
      x = (double(i)+0.5)/double(histo_bins)*(range_max-range_min) + range_min;
      y = histo[i];
    }
  }
  
  delete[] histo;
  return x;
};
void Util_CreateFakeLogLogHistogram(double** x, double** y, double** w, const double* data, const long samples, const double range_min, const double range_max, const int histo_bins)
{
  // create and fill histogram
  long* histo = new long[histo_bins];
  memset(histo, 0, histo_bins*sizeof(long));
  for(long t=0; t<samples; t++)
    histo[int(discretize(log(data[t]),log(range_min),log(range_max),histo_bins))] += 1;

  *x = new double[histo_bins];
  *y = new double[histo_bins];
  *w = new double[histo_bins];
  for(int i=0; i<histo_bins; i++)
  {
    // set up x
    if(histo[i]>2) (*y)[i] = log(double(histo[i]));
    else (*y)[i] = 0.;
    
    // set up x and w
    (*x)[i] = log(double(i+1)/double(histo_bins)*(range_max-range_min) + range_min);
    (*w)[i] = exp(-pow((double(i)-double(histo_bins)/2.)/(double(histo_bins)/6.),2.));
  }
  // normalize weights
  for(int i=0; i<histo_bins; i++)
    (*w)[i] = ((*w)[i])/total(*w,histo_bins);
  
  delete[] histo;
};
void Util_FreeFakeHistogram(double* x, double* y, double* w)
{
  delete[] x;
  delete[] y;
  delete[] w;
};

void Util_CoordinatedForMathematica(double*x, double*y, int length, IOSTREAMH)
{
  IOSTREAMC <<"{";
  for(int i=0; i<length; i++)
  {
    if(i>0) IOSTREAMC <<",";
    IOSTREAMC <<"{"<<x[i]<<","<<y[i]<<"}";
  }
  IOSTREAMC <<"};"<<IOSTREAMENDL;
};
