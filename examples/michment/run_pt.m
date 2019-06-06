% this script will start-up a parallel tempering job for
%  the Michaelis-Menten model 
%
%   to execute from command line on remote server:
%   > nohup matlab -r run_pt > run_pt.out 2> run_pt.err < /dev/null &
  
% path to parallel tempering scripts
addpath('~/Documents/GitHub/ptempest/core/');

% path to supplementary distributions
addpath('~/Documents/GitHub/ptempest/core/distr/');

% path to analysis tools (optional)
addpath('~/Documents/GitHub/ptempest/vis/');

% path to model-specific files (if not the current directory)
%addpath('~/googlecode/ptempest/examples/michment/');

% load configuration file
cfg = config();

% start parallel tempering
parallel_tempering(cfg);

% quit
quit;

