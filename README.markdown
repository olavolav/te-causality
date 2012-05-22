# TE-Causality

This package consists of the full source code used for the Transfer Entropy (TE) based measure of effective connectivity (called Generalized Transfer Entropy, or GTE) published in [Model-free Reconstruction of Neuronal Network Connectivity from Calcium Imaging Signals](http://arxiv.org/abs/1201.0732).

It also contains three methods of estimating TE and GTE without binning:

- te-binless-Leonenko, based on: L.F. Kozachenko and N.N. Leonenko, [Sample Estimate of the Entropy of a Random Vector](http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=ppi&paperid=797&option_lang=eng), Problems of Information Transmission, 1987
- te-binless-Kraskov, based on: A. Kraskov et al., [Estimating Mutual Information](http://pre.aps.org/abstract/PRE/v69/i6/e066138), Physical Review E, 2004
- te-binless-Frenzel, based on: S. Frenzel and B. Pompe, [Partial Mutual Information for Coupling Analysis of Multivariate Time Series](http://prl.aps.org/abstract/PRL/v99/i20/e204101), Physical Review Letters, 2007



## Dependencies

### GTE

To compile the GTE binaries based on binned estimates, all you need is a standard C++ compiler, and the following libraries:

- [GNU Scientific Library (GSL)](http://www.gnu.org/s/gsl/)
- [Boost](http://www.boost.org/)
- and Christoph Kirsts's [SimKernel](https://github.com/ChristophKirst/SimKernel) package

Please make sure that GSL and Boost are available to your C++ linker, and that the path to the SimKernel files is correctly set in the Rakefile.

### Light scattering

To simulate light scattering, we need to load the spatial positions of each node from a YAML file. You therefore need to have the [yaml-cpp](http://code.google.com/p/yaml-cpp) libraries installed.

### Binless methods

To use the binless estimators, you also need to install Marius Muja's excellent [FLANN](https://github.com/mariusmuja/flann) (Fast Library for Approximate Nearest Neighbors) package and make it available to your linker.



## Input file formats

See below for details on the file formats of the input data to the causality programs. You will need a SimKernel control file, and and some input time series, either directly loaded from disk, or in the form of spike trains. In the latter case, the te-datainit library will simulate a fluorescence signal.

Note that minimal sample files are included in the transferentropy/tests directory.

### SimKernel control file

All of the resconstruction and signal parameters are set via SimKernel. See the [SimKernel repository](https://github.com/ChristophKirst/SimKernel) for general documentation.

In principle, all parameters like `size` for the number of nodes or `SourceMarkovOrder` for the order of the assumed Markov process are set in this file and can be changed without re-compiling the GTE programs.

Example control file:

```c++
// word length
p = 2; 
SourceMarkovOrder = p;
TargetMarkovOrder = p;
StartSampleIndex = p;

bins = 3;

// input data
size = 100;
spikeindexfile = "data/indices.txt";
spiketimesfile = "data/times.txt";
tauF = 20;

samples = 1000;

FluorescenceModel = "Leogang";
DeltaCalciumOnAP = 50.;
tauCa = 1000.;
noise = 0.03;
saturation = 300;
tauF = 20;

globalbins = 2;
OverrideRescalingQ = False;
HighPassFilterQ = True;
InstantFeedbackTermQ = True;
IncludeGlobalSignalQ = True;
GlobalConditioningLevel = 0.5;

// output files
outputfile = "adjA.mx";
outputparsfile = "pars.mx";
```

This couls be a standard setting to simulate a fluorescence signal based on the given spike times, and the calculate the GTE matrix conditioned on the average fluorescence via the parameter `GlobalConditioningLevel`.

One of the main features of SimKernel are iterators. For instance, to calculate GTE for different numbers of samples, one could use code such as the following:

```c++
samplesList = {10,50,100,500,1000,5000,10000,15000};
iS = Iterator[j,{j,0,Length[samplesList]-1,1}];
samples = samplesList[[iS]];
```

Note the double brackets (Mathematica syntax), but that the array index is starting from zero.

### Spike trains

Spike times are loaded simply as linebreak-separated entries in two files of equal length. One containing spike times (SimKernel parameter `spiketimesfile`) and one listing the corresponding node indices (SimKernel parameter `spikeindexfile`).

### Fluorescence time series (or any other continuous signal)

Alternatively, you can simply load an already existing time series (in ASCII format) directly via the SimKernel parameter `inputfile`.

### Network topology

For the puropse of simulating light scattering, we need to know the spatial position of each node. We therefore load a YAML file of the following format:

```yaml
---
size: 3 # number of nodes
cons: 4 # number of connections
notes: "example miniature topology"
createdAt: "Tue 20 Sep 2011 14:55:21"
nodes:
  - id: 1
    pos: [0.516609, 0.187339]
    connectedTo: [2,3]
  - id: 2
    pos: [0.933885, 0.0615908]
    connectedTo: [1]
  - id: 3
    pos: [0.519384, 0.110653]
    connectedTo: [2]
```


## Copyright

All of the files can be copied, modified, used for commercial or non-commercial purpose, as long as you keep the following copyright message in the source files:

(ADD SOMETHING HERE LATER!)
