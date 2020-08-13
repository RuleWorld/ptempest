# Tools to simplify the process of starting the ptempest algorithm and visualizing the results. 

Two main functions exists within the folder: 
1. run_ptempest
2. plotting_tool

### The tools require Matlab 2020a 

## run_ptempest 
The goal of the run_ptempest function is to simplify the process of utilizing ptempest to fit models, built in BioNetGen, to data. 

## plotting_tool 
The plotting_tool function is a tool to quickly visualize the quality of the fit after running the run_ptempest program. To use the tool, the function must be instantiated where the minimum input must include the path of the output from the run_ptempest script. For example: 
~~~
x = plotting_tool(<path_to_run_ptempeset_results>) 
~~~

Instianting the tool with the minimum input, the following is assumed: 
1. The number of samples chosen is 100 (n_samples) 
2. The iteration in which ptempest converged is assumed to be at the start of the simulation, therefore (converged_start = 1) 
3. The last iteration one wishes to sample is assumed to be the last avaiable swap recorded in ptempest (end_sampling). 
