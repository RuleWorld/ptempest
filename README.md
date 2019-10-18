# PTempEst - Parallel Tempering for Estimation

*The 'P' is silent.*
Justin S. Hogg (but the rest of us call it 'pee-tempest'...)
 
<pre>

directory structure
-------------------
ptempest/             : root PTempEst directory.
ptempest/core/        : core m-files for parallel tempering.
ptempest/core/distr/  : supplementary probability distribution functions
ptempest/vis/         : scripts for visualizing posterior parameter and
                          trajectory distributions.
ptempest/examples/    : example model configurations. 
ptempest/mex_cvode/   : C language framework for implementing dynamical system
                          models  using the SUNDIALS CVODE integrator with a
                          Matlab MEX interface.
</pre>

[Instructions on setting up CVODE libraries to link to BioNetGen-generated MEX code.](https://docs.google.com/document/d/1sVKyIkFPAhjOpuexq_NEV8lUGL7pRuX_gF7BmThsUls/edit?usp=sharing)
