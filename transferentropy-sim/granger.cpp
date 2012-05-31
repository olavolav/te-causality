// calculate the Granger causality between a numer of time series
// (new version usting the te-datainit library)
// created by olav, 6 Sep 2011 17:03:14 CEST

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <ctime>
#include <cstring>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_vector_double.h>

#include "../olav.h"
#include "../../../Sonstiges/SimKernel/sim_main.h"
#include "../te-datainit.h"

#undef NORMALIZE_GRANGER_CAUSALITY_ESTIMATE

// #define GSL_RANDOM_NUMBER_GENERATOR gsl_rng_default
#define GSL_RANDOM_NUMBER_GENERATOR gsl_rng_ranlxs2

#define REPORTS 25
#undef SHOW_DETAILED_PROGRESS

#define SEPARATED_OUTPUT

using namespace std;

typedef unsigned char rawdata;

time_t start, middle, end, now;

class Kernel;

int main(int argc, char* argv[])
{
  SimControl<Kernel> simc;
  time(&start);
  int ret = simc.simulate(argc, argv);
  return 0;
};

class Kernel
{
public:
  unsigned long iteration;
  unsigned int size;
  // unsigned int bins
  unsigned int globalbins;
  // unsigned long mag der SimKernel irgendwie nicht?
  long samples, effectivesamples;
  long StartSampleIndex, EndSampleIndex;
  bool EqualSampleNumberQ;
  long MaxSampleNumberPerBin;
  unsigned long * AvailableSamples;
  double std_noise;
  string inputfile_name;
  string outputfile_results_name;
  string outputfile_pars_name;
  string spikeindexfile_name, spiketimesfile_name;
  string FluorescenceModel;

  double input_scaling;
  double cutoff;
  double tauF;
  double tauCa;
  double fluorescence_saturation;
  double DeltaCalciumOnAP;

  // parameters for light scattering
  std::string YAMLfilename;
  double SigmaScatter;
  double AmplitudeScatter;
  
  bool OverrideRescalingQ; // if, for example the input are pure spike data (integers)
  bool HighPassFilterQ; // actually, this transforms the signal into the difference signal
  bool InstantFeedbackTermQ;

  bool IncludeGlobalSignalQ;
  bool GenerateGlobalFromFilteredDataQ;
  double GlobalConditioningLevel;
  int SourceMarkovOrder, TargetMarkovOrder, MarkovOrder;
  
  bool ContinueOnErrorQ;
  bool skip_the_rest;
  
  // bool AutoBinNumberQ;
  bool AutoConditioningLevelQ;
  
  gsl_rng* GSLrandom;

  // declare GSL workspaces
  gsl_multifit_linear_workspace * GSLworkspaceBoth;
  double residue;
  gsl_matrix* inputBoth;
  gsl_matrix* covBoth;
  gsl_vector* output;
  gsl_vector* coeffBoth;

  gsl_multifit_linear_workspace * GSLworkspaceSingle;
  gsl_matrix* inputSingle;
  gsl_matrix* covSingle;
  gsl_vector* coeffSingle;

#ifdef NORMALIZE_GRANGER_CAUSALITY_ESTIMATE
  long* shuffle_permutation;
#endif

  double **xdatadouble;
  rawdata *xglobal;
#ifndef SEPARATED_OUTPUT
  double **xresult;
#else
  double ***xresult;
#endif

  void initialize(Sim& sim)
  {
    iteration = sim.iteration();
    sim.io <<"Init: iteration "<<iteration<<", process "<< sim.process()<<Endl;
    time(&now);
    sim.io <<"time: ";
    sim.io <<"elapsed "<<sec2string(difftime(now,start));
    sim.io <<", ETA "<<ETAstring(sim.iteration()-1,sim.n_iterations(),difftime(now,start))<<Endl;
    
    // read parameters from control file
    sim.get("size",size);
    // bins = 0;
    // sim.get("AutoBinNumberQ",AutoBinNumberQ,false);
    // if(!AutoBinNumberQ) sim.get("bins",bins);
    sim.get("globalbins",globalbins,1);
    if(globalbins<1) globalbins=1;
    sim.get("samples",samples);
    sim.get("EqualSampleNumberQ",EqualSampleNumberQ,false);
    sim.get("MaxSampleNumberPerBin",MaxSampleNumberPerBin,-1);

    sim.get("noise",std_noise,-1.);
    if(std_noise<0) {
      sim.io <<"Noise hast to be positive for the NN estimator to work. Aborting."<<Endl;
      exit(1);
    }
    sim.get("appliedscaling",input_scaling,1.);
    sim.get("cutoff",cutoff,-1.);
    sim.get("tauF",tauF);
    sim.get("OverrideRescalingQ",OverrideRescalingQ,false);
    sim.get("HighPassFilterQ",HighPassFilterQ,false);
    sim.get("InstantFeedbackTermQ",InstantFeedbackTermQ,false);

    sim.get("saturation",fluorescence_saturation,-1.);
    sim.get("IncludeGlobalSignalQ",IncludeGlobalSignalQ,false);
    assert(IncludeGlobalSignalQ ^ globalbins==1);
    sim.get("GenerateGlobalFromFilteredDataQ",GenerateGlobalFromFilteredDataQ,false);
    sim.get("AutoConditioningLevelQ",AutoConditioningLevelQ,false);
    if(!AutoConditioningLevelQ) {
      sim.get("GlobalConditioningLevel",GlobalConditioningLevel,-1.);
      if (GlobalConditioningLevel>0) assert(globalbins==2); // for now, to keep it simple
    }
    else GlobalConditioningLevel = -1.;
    
    sim.get("SourceMarkovOrder",SourceMarkovOrder,1);
    assert(SourceMarkovOrder>0);
    sim.get("TargetMarkovOrder",TargetMarkovOrder,1);
    assert(TargetMarkovOrder==SourceMarkovOrder);
    MarkovOrder = SourceMarkovOrder; // we will here assume equal Markov order for ease of use

    sim.get("StartSampleIndex",StartSampleIndex,1);
    assert(StartSampleIndex>=MarkovOrder && StartSampleIndex<samples);
    sim.get("EndSampleIndex",EndSampleIndex,samples-1);
    assert(EndSampleIndex>=StartSampleIndex && EndSampleIndex<samples);
    effectivesamples = EndSampleIndex-StartSampleIndex+1;
    
    sim.get("inputfile",inputfile_name,"");
    sim.get("outputfile",outputfile_results_name);
    sim.get("outputparsfile",outputfile_pars_name);
    sim.get("spikeindexfile",spikeindexfile_name,"");
    sim.get("spiketimesfile",spiketimesfile_name,"");
    // make sure we either have a fluorescence input or both spike inputs
    if(!((inputfile_name!="") ^ ((spikeindexfile_name!="")&&(spiketimesfile_name!="")))) {
      sim.io <<"Error: Based on the parameters, it is not clear where your data should come from."<<Endl;
      exit(1);
    }
    sim.get("FluorescenceModel",FluorescenceModel,"");
    sim.get("DeltaCalciumOnAP",DeltaCalciumOnAP,50);
    sim.get("tauCa",tauCa,1000);

    // parameters for light scattering
    sim.get("YAMLfile",YAMLfilename,"");
    sim.get("SigmaScatter",SigmaScatter,-1.);
    sim.get("AmplitudeScatter",AmplitudeScatter,-1.);
  
    sim.get("ContinueOnErrorQ",ContinueOnErrorQ,false);
    
    // initialize random number generator
    gsl_rng_env_setup();
    GSLrandom = gsl_rng_alloc(GSL_RANDOM_NUMBER_GENERATOR);
    gsl_rng_set(GSLrandom, 1234);
    
    // declare GSL workspaces
    GSLworkspaceBoth = NULL;
    inputBoth = NULL;
    covBoth = NULL;
    output = NULL;
    coeffBoth = NULL;

    GSLworkspaceSingle = NULL;
    inputSingle = NULL;
    covSingle = NULL;
    coeffSingle = NULL;

#ifdef NORMALIZE_GRANGER_CAUSALITY_ESTIMATE
    shuffle_permutation = NULL;
#endif

    AvailableSamples = NULL;
    xdatadouble = NULL;
    xglobal = NULL;
    xresult = NULL;
  };

  void execute(Sim& sim)
  {
    sim.io <<"------ granger-sim:v2 ------ olav, 06 Sep 2011 ------"<<Endl;
    // time_t start, middle, end;

    sim.io <<"output file: "<<outputfile_results_name<<Endl;
    // Gespeichert wird später - hier nur Test, ob das Zielverzeichnis existiert
    write_parameters();
    skip_the_rest = false;

    sim.io <<"allocating memory..."<<Endl;
    try {
#ifndef SEPARATED_OUTPUT
      xresult = new double*[size];
#else
      xresult = new double**[size];
#endif
      for(int i=0; i<size; i++) {
#ifndef SEPARATED_OUTPUT
        xresult[i] = new double[size];
        memset(xresult[i], 0, size*sizeof(double));
#else
        xresult[i] = new double*[size];
        for(int i2=0; i2<size; i2++)
        {
          xresult[i][i2] = new double[globalbins];
          memset(xresult[i][i2], 0, globalbins*sizeof(double));
        }
#endif
      }
    
      AvailableSamples = new unsigned long[globalbins];
            
      // hack of medium ugliness to make it work without global signal
      if(globalbins<=1) {
        // sim.io <<"debug: xglobal hack."<<Endl;
        xglobal = new rawdata[samples];
        memset(xglobal, 0, samples*sizeof(rawdata));
        AvailableSamples[0] = EndSampleIndex-StartSampleIndex+1;
      }
      sim.io <<" -> done."<<Endl;     
      
      // double** xdatadouble = NULL;
      if(inputfile_name!="")
        sim.io <<"input file: \""<<inputfile_name<<"\""<<Endl;
      else {
        sim.io <<"input files:"<<Endl;
        sim.io <<"- spike indices: \""<<spikeindexfile_name<<"\""<<Endl;
        sim.io <<"- spike times: \""<<spiketimesfile_name<<"\""<<Endl;
      }
      
      if(inputfile_name=="") {
        sim.io <<"loading data and generating time series from spike data..."<<Endl;
        xdatadouble = generate_time_series_from_spike_data(spiketimesfile_name, spikeindexfile_name, size, int(round(tauF)), samples, FluorescenceModel, std_noise, fluorescence_saturation, cutoff, DeltaCalciumOnAP, tauCa, GSLrandom, sim);
      }
      else {
        sim.io <<"loading data from binary file..."<<Endl;
        xdatadouble = load_time_series_from_file(inputfile_name, size, samples, input_scaling, OverrideRescalingQ, std_noise, fluorescence_saturation, cutoff, GSLrandom, sim);
      }
      sim.io <<" -> done."<<Endl;
      
      if(AmplitudeScatter>0.) {
        sim.io <<"simulating light scattering..."<<Endl;
        apply_light_scattering_to_time_series(xdatadouble, size, samples, YAMLfilename, SigmaScatter, AmplitudeScatter, sim);
        sim.io <<" -> done."<<Endl;
      }

      // sim.io <<"histogram of averaged signal:"<<Endl;
      // double* xmean = generate_mean_time_series(xdatadouble,size,samples);
      // PlotLogHistogramInASCII(xmean,samples,smallest(xmean,samples),largest(xmean,samples),"<fluoro>","#",sim);
      // free_time_series_memory(xmean);

      if(AutoConditioningLevelQ) {
        sim.io <<"guessing optimal conditioning level..."<<Endl;
        GlobalConditioningLevel = Magic_GuessConditioningLevel(xdatadouble,size,samples,sim);
        sim.io <<" -> conditioning level is: "<<GlobalConditioningLevel<<Endl;
        sim.io <<" -> done."<<Endl;
      }      
      
      if((globalbins>1)&&(!GenerateGlobalFromFilteredDataQ)) {
        sim.io <<"generating discretized global signal..."<<Endl;
        xglobal = generate_discretized_global_time_series(xdatadouble, size, samples, globalbins, GlobalConditioningLevel, AvailableSamples, StartSampleIndex, EndSampleIndex, EqualSampleNumberQ, MaxSampleNumberPerBin, sim);
        sim.io <<" -> done."<<Endl;
      }
      
      if(HighPassFilterQ) {
        sim.io <<"applying high-pass filter to time series..."<<Endl;
        apply_high_pass_filter_to_time_series(xdatadouble, size, samples);
        sim.io <<" -> done."<<Endl;
      }
      
      if((globalbins>1)&&(GenerateGlobalFromFilteredDataQ)) {
        sim.io <<"generating discretized global signal..."<<Endl;
        xglobal = generate_discretized_global_time_series(xdatadouble, size, samples, globalbins, GlobalConditioningLevel, AvailableSamples, StartSampleIndex, EndSampleIndex, EqualSampleNumberQ, MaxSampleNumberPerBin, sim);
        sim.io <<" -> done."<<Endl;
      }

    // allocate GSL workspaces
    GSLworkspaceBoth = gsl_multifit_linear_alloc(samples-MarkovOrder,2*MarkovOrder+1);
    // the coulmns in input are: (nr. of columns)
    // target (MarkovOrder), source (MarkovOrder), constant term (1)
    inputBoth = gsl_matrix_alloc(samples-MarkovOrder,2*MarkovOrder+1);
    // because the last column of this matrix never changes, we can set it here:
    gsl_matrix_set_all(inputBoth,1.0);
    covBoth = gsl_matrix_alloc(2*MarkovOrder+1,2*MarkovOrder+1);
    // output is the only GSL consturuct that we use for both cases ("both" and "single")
    output = gsl_vector_alloc(samples-MarkovOrder);
    coeffBoth = gsl_vector_alloc(2*MarkovOrder+1);

    GSLworkspaceSingle = gsl_multifit_linear_alloc(samples-MarkovOrder,1*MarkovOrder+1);
    inputSingle = gsl_matrix_alloc(samples-MarkovOrder,1*MarkovOrder+1);
    gsl_matrix_set_all(inputSingle,1.0);
    covSingle = gsl_matrix_alloc(1*MarkovOrder+1,1*MarkovOrder+1);
    coeffSingle = gsl_vector_alloc(1*MarkovOrder+1);


#ifdef NORMALIZE_GRANGER_CAUSALITY_ESTIMATE
      // shuffle_permutation = generate_random_permutation(samples,globalbins,AvailableSamples,StartSampleIndex,EndSampleIndex,xglobal);
      shuffle_permutation = generate_random_geometric_permutation(samples,globalbins,xglobal,5,GSLrandom);
#endif
    }
    catch(...) {
      sim.io <<"Error: could not reserve enough memory!"<<Endl;
      if(!ContinueOnErrorQ) exit(1);
      else {
        sim.io <<"Error handling: ContinueOnErrorQ flag set, proceeding..."<<Endl;
        skip_the_rest = true;
      }
    }
    // sim.io <<" -> done."<<Endl;
      
    if (!skip_the_rest) {   
      // ----------------------------------------- main loop: start -----------------------------------------
      sim.io <<"set-up: "<<size<<" nodes, ";
      sim.io <<EndSampleIndex-StartSampleIndex+1<<" out of "<<samples<<" samples, ";
      // sim.io <<bins<<" bins, "
      sim.io <<globalbins<<" globalbins"<<Endl;
      sim.io <<"set-up: Markov order of source/target/conditioning: "<<SourceMarkovOrder<<"/"<<TargetMarkovOrder<<"/1"<<Endl;
#ifdef SEPARATED_OUTPUT
      sim.io <<"set-up: separated output (globalbin)"<<Endl;
#endif
  
      time(&start);
      sim.io <<"start: "<<ctime(&start)<<Endl;
#ifdef SHOW_DETAILED_PROGRESS
      sim.io <<"running ";
#else
      sim.io <<"running..."<<Endl;
      bool status_already_displayed = false;
#endif

      for(int ii=0; ii<size; ii++)
      {
#ifdef SHOW_DETAILED_PROGRESS
        status(ii,REPORTS,size);
#else
        time(&middle);
        if ((!status_already_displayed)&&((ii>=size/3)||((middle-start>30.)&&(ii>0)))) { 
          sim.io <<" (after "<<ii<<" nodes: elapsed "<<sec2string(difftime(middle,start)) \
            <<", ETA "<<ETAstring(ii,size,difftime(middle,start))<<")"<<Endl;
          status_already_displayed = true;
        }
#endif
        for(int jj=0; jj<size; jj++)
        {
          if (ii != jj) {
#ifndef SEPARATED_OUTPUT
            xresult[jj][ii] = GrangerCausality(xdatadouble[ii], xdatadouble[jj]);
#else
            GrangerCausalitySeparated(xdatadouble[ii], xdatadouble[jj], ii, jj);
#endif
          }
          // else xresult[ii][jj] = 0.0;
        }
      }
#ifndef SHOW_DETAILED_PROGRESS
    sim.io <<" -> done."<<Endl;
#endif
    }
    // ----------------------------------------- main loop: end -----------------------------------------
    time(&end);
    sim.io <<"end: "<<ctime(&end)<<Endl;
    sim.io <<"runtime: "<<sec2string(difftime(end,start))<<Endl;

    // cout <<"TE terms: "<<terms_sum<<", of those zero: "<<terms_zero<<" ("<<int(double(terms_zero)*100/terms_sum)<<"%)"<<endl;
  };
  
  void finalize(Sim& sim) {
    if(!skip_the_rest) {
#ifdef SEPARATED_OUTPUT
      write_multidim_result(xresult,globalbins,size,outputfile_results_name,sim);
#else
      write_result(xresult,size,outputfile_results_name,sim);
#endif
      write_parameters();

      // free allocated memory
      gsl_rng_free(GSLrandom);  
#ifdef NORMALIZE_GRANGER_CAUSALITY_ESTIMATE
      if (shuffle_permutation != NULL ) delete[] shuffle_permutation;
#endif
      
    // free memory
    gsl_multifit_linear_free(GSLworkspaceBoth);
    gsl_matrix_free(inputBoth);
    gsl_matrix_free(covBoth);
    gsl_vector_free(output);
    gsl_vector_free(coeffBoth);
    gsl_multifit_linear_free(GSLworkspaceSingle);
    gsl_matrix_free(inputSingle);
    gsl_matrix_free(covSingle);
    gsl_vector_free(coeffSingle);
    }

    try {
#ifdef SEPARATED_OUTPUT
      for (int x=0; x<size; x++) {
        for (int x2=0; x2<size; x2++)
          delete[] xresult[x][x2];
        delete[] xresult[x];
      }
#endif
      delete[] xresult;
    
      if (AvailableSamples != NULL) delete[] AvailableSamples;

      if(xdatadouble != NULL) free_time_series_memory(xdatadouble,size);
      if(xglobal != NULL) free_time_series_memory(xglobal);
    }
    catch(...) {};
    
    sim.io <<"End of Kernel (iteration="<<(sim.iteration())<<")"<<Endl;
  };


  void write_parameters() {
    char* name = new char[outputfile_pars_name.length()+1];
    strcpy(name,outputfile_pars_name.c_str());
    ofstream fileout1(name);
    delete[] name;
    if (fileout1 == NULL) {
      cerr <<endl<<"error: cannot open output file!"<<endl;
      exit(1);
    }

    fileout1.precision(6);
    fileout1 <<"{";
    fileout1 <<"executable->granger";
#ifdef NORMALIZE_GRANGER_CAUSALITY_ESTIMATE
    fileout1 <<", NormalizationViaShufflingQ->True";
#else
    fileout1 <<", NormalizationViaShufflingQ->False";
#endif
    fileout1 <<", iteration->"<<iteration;
    time(&end);
    fileout1 <<", ExecutionTime->\""<<sec2string(difftime(end,start))<<"\"";
    
    fileout1 <<", size->"<<size;
    // fileout1 <<", bins->"<<bins;
    fileout1 <<", globalbins->"<<globalbins;
    fileout1 <<", appliedscaling->"<<input_scaling;
    fileout1 <<", samples->"<<samples;
    fileout1 <<", StartSampleIndex->"<<StartSampleIndex;
    fileout1 <<", EndSampleIndex->"<<EndSampleIndex;
    fileout1 <<", EqualSampleNumberQ->"<<bool2textMX(EqualSampleNumberQ);
    fileout1 <<", MaxSampleNumberPerBin->"<<MaxSampleNumberPerBin;
    fileout1 <<", AvailableSamples->{";
    for (int i=0; i<globalbins; i++) {
      if (i>0) fileout1 <<",";
      if (AvailableSamples == NULL) fileout1 <<"?";
      else fileout1 <<AvailableSamples[i];
    }
    fileout1 <<"}";

    fileout1 <<", cutoff->"<<cutoff;
    fileout1 <<", noise->"<<std_noise;
    fileout1 <<", tauF->"<<tauF;
    fileout1 <<", tauCa->"<<tauCa;
    fileout1 <<", DeltaCalciumOnAP->"<<DeltaCalciumOnAP;
    fileout1 <<", OverrideRescalingQ->"<<bool2textMX(OverrideRescalingQ);
    fileout1 <<", HighPassFilterQ->"<<bool2textMX(HighPassFilterQ);
    fileout1 <<", InstantFeedbackTermQ->"<<bool2textMX(InstantFeedbackTermQ);
    fileout1 <<", ContinueOnErrorQ->"<<bool2textMX(ContinueOnErrorQ);
    fileout1 <<", saturation->"<<fluorescence_saturation;
    fileout1 <<", IncludeGlobalSignalQ->"<<bool2textMX(IncludeGlobalSignalQ);
    fileout1 <<", GenerateGlobalFromFilteredDataQ->"<<bool2textMX(GenerateGlobalFromFilteredDataQ);
    fileout1 <<", GlobalConditioningLevel->"<<GlobalConditioningLevel;
    fileout1 <<", TargetMarkovOrder->"<<TargetMarkovOrder;
    fileout1 <<", SourceMarkovOrder->"<<SourceMarkovOrder;
    
    // fileout1 <<", AutoBinNumberQ->"<<bool2textMX(AutoBinNumberQ);
    fileout1 <<", AutoConditioningLevelQ->"<<bool2textMX(AutoConditioningLevelQ);
    
    
    fileout1 <<", inputfile->\""<<inputfile_name<<"\"";
    fileout1 <<", outputfile->\""<<outputfile_results_name<<"\"";
    fileout1 <<", spikeindexfile->\""<<spikeindexfile_name<<"\"";
    fileout1 <<", spiketimesfile->\""<<spiketimesfile_name<<"\"";
    fileout1 <<", FluorescenceModel->\""<<FluorescenceModel<<"\"";
    // parameters for light scattering
    fileout1 <<", YAMLfile->\""<<YAMLfilename<<"\"";
    fileout1 <<", SigmaScatter->"<<SigmaScatter;
    fileout1 <<", AmplitudeScatter->"<<AmplitudeScatter;
    fileout1 <<"}"<<endl;
    
    fileout1.close();
  };
  

#ifdef SEPARATED_OUTPUT
  void GrangerCausalitySeparated(double *arrayI, double *arrayJ, int I, int J)
#else
  double GrangerCausality(double *arrayI, double *arrayJ)
#endif
  {
    // examining GC J -> I
    long t_sample;
    double result;
    long const JShift = (long const)InstantFeedbackTermQ;
    
    for(rawdata g=0; g<globalbins; g++) {
      if(AvailableSamples[g] > 1) {
        t_sample = 0;
        for(long t=StartSampleIndex; t<=EndSampleIndex; t++) {
          if(xglobal[t] == g) {
            // 1.) set up source and target data of connection
            for(long tt2=0; tt2<MarkovOrder; tt2++) {
              gsl_matrix_set(inputSingle,t_sample,tt2,arrayI[t-1-tt2]);
              gsl_matrix_set(inputBoth,t_sample,tt2,arrayJ[t-1-tt2 + JShift]);
              gsl_matrix_set(inputBoth,t_sample,MarkovOrder+tt2,arrayI[t-1-tt2]);
            }
            gsl_vector_set(output,t_sample,arrayI[t]);
            
            t_sample++;
          }
        }
        assert(t_sample==AvailableSamples[g]);

        // 2.) calculate variance of residue of single time series GC
        gsl_multifit_linear(inputSingle,output,coeffSingle,covSingle,&residue,GSLworkspaceSingle);
        // test_residue = residue;
        
        // 3.) calculate variance of residue of both time series GC
        gsl_multifit_linear(inputBoth,output,coeffBoth,covBoth,&residue,GSLworkspaceBoth);

        // test:
        // xresult[ii][jj] = residue/test_residue;

        // save those coefficients that depend on the source
        // covTraceBoth = covTraceSingle = 0.0;
        // for (int kk=0; kk<MarkovOrder; kk++)
        // {
          // xresult[ii][jj][kk] = gsl_vector_get(coeffBoth,MarkovOrder+kk);
          // covTraceBoth += gsl_matrix_get(covBoth,MarkovOrder+kk,MarkovOrder+kk);
          // covTraceSingle += gsl_matrix_get(covSingle,MarkovOrder+kk,MarkovOrder+kk);
        // }
        // xresult[ii][jj][MarkovOrder] = residue;
        // xresult[ii][jj] = gsl_vector_get(coeffBoth,MarkovOrder);

        // save reduction in variance as xresult
        result = \
          log(sqrt(gsl_matrix_get(covSingle,0,0)+gsl_matrix_get(covSingle,1,1)) / \
          sqrt(gsl_matrix_get(covBoth,0,0)+gsl_matrix_get(covBoth,1,1)));
        // result = \
        //   sqrt(gsl_matrix_get(covBoth,0,0)+gsl_matrix_get(covBoth,1,1)) / \
        //   sqrt(gsl_matrix_get(covSingle,0,0)+gsl_matrix_get(covSingle,1,1)); // ask Demian about this!
          
#ifdef SEPARATED_OUTPUT
        xresult[I][J][g] = result;
#endif        
      }
    }
#ifndef SEPARATED_OUTPUT
    return result;
#endif
  };
  
  
};

