// created Do 28 Apr 2011 18:13:44 CEST by olav

#include "te-datainit.h"

#define FMODEL_SPIKECOUNT 1
#define FMODEL_HOWMANYAREACTIVE 2
#define FMODEL_LEOGANG 3
#define FMODEL_ERROR -1

double** load_time_series_from_binary_file(std::string inputfile_name, unsigned int size, long samples, double input_scaling, bool OverrideRescalingQ, double std_noise, double fluorescence_saturation, double cutoff)
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
  	std::cout <<std::endl<<"error: cannot find input file!"<<std::endl;
  	exit(1);
  }

	// test file length
	binaryfile.seekg(0,std::ios::end);
	if(long(binaryfile.tellg()) != size*samples)
	{
  	std::cout <<std::endl<<"error: file length of input does not match given parameters!"<<std::endl;
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
  //  std::cout <<"DEBUG: maxsamples = "<<maxsamples<<std::endl;
  //  
  //  if ((MaxSampleNumberPerBin>maxsamples) && !EqualSampleNumberQ)
  //    maxsamples = MaxSampleNumberPerBin;
  //  if ((MaxSampleNumberPerBin<maxsamples)&&(MaxSampleNumberPerBin>0))
  //    maxsamples = MaxSampleNumberPerBin;
  //  std::cout <<"DEBUG: cut to maxsamples = "<<maxsamples<<std::endl;
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
	
  delete[] in_from_file_array;
  // delete[] xglobaltemp;
	delete[] tempdoublearray;
  // delete[] tempdoublearraycopy;
	
  return xresult;
};


rawdata* generate_discretized_global_time_series(double** time_series, unsigned int size, long samples, unsigned int globalbins, double GlobalConditioningLevel)
{
	rawdata* xglobal = new rawdata[samples];
	memset(xglobal, 0, samples*sizeof(rawdata));
	double* xglobaltemp = new double[samples];
	memset(xglobaltemp, 0, samples*sizeof(double));

	for (long t=0; t<samples; t++)
  {
    for (unsigned int ii=0; ii<size; ii++)
  		xglobaltemp[t] += time_series[ii][t];
		xglobaltemp[t] /= size;
  }
  	
	// EVIL SAALBACH HACK FOR TIME CODE GLOBAL SIGNAL: -------------------------------------------- !!!!!!!!!!
	// xglobaltemp[0] = 0.;
	// for (unsigned long t=0; t<samples; t++)
	// 	xglobaltemp[t] = double(int(t)%int(60*24/tauF));
	// 	// xglobaltemp[t] = double(mod(t,60*24/tauF));
	// std::cout <<"DEBUG OF EVIL TIME CODE HACK: last globaltemp value = "<<xglobaltemp[samples-1]<<std::endl;
	
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
		std::cout <<"global conditioning with level "<<GlobalConditioningLevel<<": "<<(100.*below)/samples<<"% are below threshold... "<<std::endl;
	}
	else discretize(xglobaltemp,xglobal,samples,globalbins);
	
	// determine available samples per globalbin for TE normalization later
  // memset(AvailableSamples, 0, globalbins*sizeof(unsigned long));
  // for (unsigned long t=StartSampleIndex; t<EndSampleIndex; t++)
  //  AvailableSamples[xglobal[t]]++;


  // if (EqualSampleNumberQ || (MaxSampleNumberPerBin>0))
  // {
  //  unsigned long maxsamples = ULONG_MAX;
  //  for (rawdata g=0; g<globalbins; g++)
  //    if (AvailableSamples[g]<maxsamples) maxsamples = AvailableSamples[g];
  //  std::cout <<"DEBUG: maxsamples = "<<maxsamples<<std::endl;
  //  
  //  if ((MaxSampleNumberPerBin>maxsamples) && !EqualSampleNumberQ)
  //    maxsamples = MaxSampleNumberPerBin;
  //  if ((MaxSampleNumberPerBin<maxsamples)&&(MaxSampleNumberPerBin>0))
  //    maxsamples = MaxSampleNumberPerBin;
  //  std::cout <<"DEBUG: cut to maxsamples = "<<maxsamples<<std::endl;
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


  delete[] xglobaltemp;
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
// void discretize(double** in, rawdata** out, unsigned int size, long nr_samples, unsigned int nr_bins)
// {
//   for(unsigned int ii=0; ii<size; ii++)
//     discretize(in[ii], out[ii], nr_samples, nr_bins);
// };
void discretize(double* in, rawdata* out, long nr_samples, unsigned int nr_bins)
{
  discretize(in,out,smallest(in,nr_samples),largest(in,nr_samples),nr_samples,nr_bins);
};
void discretize(double* in, rawdata* out, double min, double max, long nr_samples, unsigned int nr_bins)
{
	// correct discretization according to 'te-test.nb'
	double xstepsize = (max-min)/nr_bins;
	// std::cout <<"max = "<<max<<std::endl;
	// std::cout <<"min = "<<min<<std::endl;
	// std::cout <<"bins here = "<<nr_bins<<std::endl;
	// std::cout <<"stepsize = "<<xstepsize<<std::endl;

	int xint;
	for (unsigned long t=0; t<nr_samples; t++)
	{
		assert(in[t]<=max);
		assert(in[t]>=min);
		if (in[t]>=max) xint = nr_bins-1;
		else
		{
			if (in[t]<=min) xint = 0;
			else xint = (int)((in[t]-min)/xstepsize); // if(!(xint<nr_bins)) std::cout<<"!"<<xint<<","<<in[t]<<std::endl; }
		}
		if (xint >= nr_bins) xint = nr_bins-1; // need to have this for some silly numerical reason...
		
		assert(xint>=0);

		out[t] = (rawdata)xint;
		assert(out[t]<nr_bins);
	}
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

double** generate_time_series_from_spike_data(std::string inputfile_spiketimes, std::string inputfile_spikeindices, unsigned int size, unsigned int tauImg, long samples, std::string fluorescence_model, double std_noise, double fluorescence_saturation, double cutoff)
{
	std::cout <<"generating time series from spike data..."<<std::endl;
	
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
	const long nr_spikes = binaryfileI.tellg();
	std::cout <<"number of spikes in index file: "<<nr_spikes<<std::endl;
	binaryfileI.seekg(0,std::ios::beg);
	int* xindex = new int[nr_spikes];
	double* xtimes = new double[nr_spikes];
	
	// load spike data
	binaryfileI.read((char*)xindex, nr_spikes*sizeof(int));
	binaryfileT.read(reinterpret_cast<char*>(xtimes), nr_spikes*sizeof(double));
	
  for(long t=0; t<20; t++)
    std::cout <<"xindex = "<<xindex[t]<<", xtimes = "<<xtimes[t]<<std::endl;

  // generate switch key for the fluorescence model
  int fluorescence_model_key = FMODEL_ERROR;
  if (fluorescence_model == "SpikeCount") fluorescence_model_key = FMODEL_SPIKECOUNT;
  if (fluorescence_model == "HowManyAreActive") fluorescence_model_key = FMODEL_HOWMANYAREACTIVE;
  if (fluorescence_model == "Leogang") fluorescence_model_key = FMODEL_LEOGANG;
  
	// generate fluorescence data
  // const int int_tauF = (int)round(tauF); // in ms
  unsigned long startindex = 1;
  unsigned long endindex = 0; // therefore, we miss the first spike here
  long dataindex = 0;
  unsigned long tinybit_spikenumber;
  // unsigned long ttExactMS = 0;
  for (unsigned long ttExactMS=0; ttExactMS<tauImg*samples; ttExactMS+=tauImg)
  {
    // determine starting and ending spike index of current frame
    while ((endindex+1<=nr_spikes)&&(xtimes[endindex+1]<=ttExactMS+tauImg))
    {
      // std::cout <<"     +1 spike!"<<std::endl;
      endindex++;
    }
    tinybit_spikenumber = std::max((int)(endindex-startindex+1),0);
    // std::cout <<"ttExactMS = "<<ttExactMS<<", startindex = "<<startindex<<", endindex = "<<endindex;
    // std::cout <<", tinybit_spikenumber = "<<tinybit_spikenumber<<std::endl;

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
          xresult[ii][dataindex] = (1.-double(tauImg)/double(1000))*xresult[ii][std::max(dataindex-1,long(0))] + \
            50.*double(count(xindex,startindex,endindex,ii));
          break;
          
        default:
          std::cout <<"error in generate_time_series_from_spike_data: invalid fluorescence model";
          exit(1);
      }
    }
    if(startindex <= endindex)
      startindex = 1 + endindex;
    dataindex++;
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

double smallest(double* array, long length)
{
	double min = array[0];
	for (long i=1; i<length; i++)
		if(array[i]<min) min = array[i];

	return min;
};
double largest(double* array, long length)
{
	double max = array[0];
	for (long i=1; i<length; i++)
		if(array[i]>max) max = array[i];

	return max;
};
double smallest(rawdata* array, long length)
{
	rawdata min = array[0];
	for (long i=1; i<length; i++)
		if(array[i]<min) min = array[i];

	return min;
};
double largest(rawdata* array, long length)
{
	rawdata max = array[0];
	for (long i=1; i<length; i++)
		if(array[i]>max) max = array[i];

	return max;
};

void free_time_series_memory(double** xresult, unsigned int size)
{
  for(unsigned int ii=0; ii<size; ii++)
    free_time_series_memory(xresult[ii]);
  delete[] xresult;
}
void free_time_series_memory(double* xresult)
{
  delete[] xresult;
};
void free_time_series_memory(rawdata** xresult, unsigned int size)
{
  for(unsigned int ii=0; ii<size; ii++)
    free_time_series_memory(xresult[ii]);
  delete[] xresult;
}
void free_time_series_memory(rawdata* xresult)
{
  delete[] xresult;
};
