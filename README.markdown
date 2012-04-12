# TE-Causality

This package consists of the full source code used for the Transfer Entropy (TE) based measure of effective connectivity (called Generalized Transfer Entropy, or GTE) published in [Model-free reconstruction of neuronal network connectivity from calcium imaging signals](http://arxiv.org/abs/1201.0732).

It also contains three methods of estimating TE and GTE without binning, based on the following publications:

- te-binless-Leonenko, based on [Sample Estimate of the Entropy of a Random Vector](http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=ppi&paperid=797&option_lang=eng)
- te-binless-Kraskov, based on [Estimating mutual information](http://pre.aps.org/abstract/PRE/v69/i6/e066138)
- te-binless-Frenzel, based on [Partial Mutual Information for Coupling Analysis of Multivariate Time Series](http://prl.aps.org/abstract/PRL/v99/i20/e204101)



## Dependencies

### GTE

To compile the GTE binaries, all you need is a standard C++ compiler, and the following libraries:

- [GNU Scientific Library (GSL)](http://www.gnu.org/s/gsl/)
- [Boost](http://www.boost.org/)
- and Christoph Kirsts's [SimKernel](https://github.com/ChristophKirst/SimKernel) package

Please make sure that GSL and Boost are available to your C++ linker, and that the path to the SimKernel files is correctly set in the Rakefile.

### Binless methods

To use the binless estimators, you also need to install Marius Muja's excellent [FLANN](https://github.com/mariusmuja/flann) (Fast Library for Approximate Nearest Neighbors) package and make it available to your linker.



## Copyright

All of the files can be copied, modified, used for commercial or non-commercial purpose, as long as you keep the following copyright message in the source files:

(ADD SOMETHING HERE LATER!)
