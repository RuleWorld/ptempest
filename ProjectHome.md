**ptempest**, _Parallel Tempering for Estimation_ (silent p), is a collection of Matlab scripts for Bayesian parameter estimation of dynamical system models. PTempEst samples parameters using the parallel tempering algorithm (a.k.a. replica exchange), a Monte Carlo method in the Metropolis-Hastings family. PTempEst supports a variety of parameter priors, including Gaussian and Laplace priors for L1 (Lasso) and L2 (Ridge) regularization.

ptempest is distributed with a few scripts for visualizing posterior parameter and trajectory distributions.

To use ptempest, the user must provide a handle to a model simulator function, create a data structure with experimental measurements, and write/edit a configuration file.

ptempest is able to utilize the parallel toolbox features in Matlab.

Project Goals for the future:
  1. Simplify the configuration process.
  1. Integrate ptempest with [BioNetGen](http://code.google.com/p/bionetgen) (rule-based modeling platform).
  1. Implement alternative version that does not depend on Matlab (Python/NumPy or pure C++)