function [cfg] = config()
%CONFIG configure negative feedback oscillatior and parallel tempering options
%
%   [cfg] = config()
%


% Initialize default configuration options 
cfg = init_config_defaults();


%% Parameter definitions (see BNGL model for parameter descriptions) [REQUIRED]
% Available prior distributions:
%   'point'           args: 'value'
%   'uniform'         args: 'min', 'max'
%   'normal',         args: 'mu', 'sigma'
%   'lognormal',      args: 'mu', 'sigma'
%   'laplace',        args: 'mu', 'b'
%   'boundedlaplace', args: 'mu', 'b', 'min', 'max'
%   'beta',           args: 'alpha', 'beta'
%   'exponential',    args: 'mu'
scale = 1;  % global weight parameter for prior strength (or "width"). scale>0.
cfg.param_defs = { ...
  struct('name','S',       'prior','point', 'value', 1, 'units','(none)' ), ...
  struct('name','YT_0',    'prior','point', 'value', 1, 'units','uM' ), ...% 1 uM
  struct('name','RT_0',    'prior','point', 'value', 1, 'units','uM' ), ...% 1 uM
  struct('name','log_k1',  'prior','laplace', 'mu', log(0.1), 'b',scale*log(2), 'units','uM/s'  ), ...% 1 uM/s
  struct('name','log_k2',  'prior','laplace', 'mu', log(0.1), 'b',scale*log(2), 'units','/s'    ), ...% 0.01 /s
  struct('name','log_k2p', 'prior','laplace', 'mu', log(0.1), 'b',scale*log(2), 'units','/uM/s' ), ...% 10 /uM/s
  struct('name','log_k3',  'prior','laplace', 'mu', log(0.1), 'b',scale*log(2), 'units','/s'    ), ...% 0.1 /uM/s
  struct('name','log_Km3', 'prior','laplace', 'mu', log(0.1), 'b',scale*log(2), 'units','uM'    ), ...% 0.01 uM
  struct('name','log_k4',  'prior','laplace', 'mu', log(0.1), 'b',scale*log(2), 'units','uM/s'  ), ...% 0.2 /uM/s
  struct('name','log_Km4', 'prior','laplace', 'mu', log(0.1), 'b',scale*log(2), 'units','uM'    ), ...% 0.01 uM
  struct('name','log_k5',  'prior','laplace', 'mu', log(0.1), 'b',scale*log(2), 'units','/s'    ), ...% 0.1 /uM/s
  struct('name','log_Km5', 'prior','laplace', 'mu', log(0.1), 'b',scale*log(2), 'units','uM'    ), ...% 0.01 uM
  struct('name','log_k6',  'prior','laplace', 'mu', log(0.1), 'b',scale*log(2), 'units','uM/s'  ), ...% 0.05 /uM/s
  struct('name','log_Km6', 'prior','laplace', 'mu', log(0.1), 'b',scale*log(2), 'units','uM'    ), ...% 0.01 uM
};
% initialize parameter distributions
cfg = init_parameter_defs( cfg.param_defs, cfg );


%% Observable definitions [REQUIRED]
%  Set display=0 to disable observable plot during visualization. Set
%  minplot/maxplot to see y-axis bounds in visualization scripts.
%  The 'norm' field should be the name of the observable that is used to
%  normalize the observable (e.g. YP is normalized by YT, the total quantity
%  of Y). Leave this field empty if normalization is not desired. NOTE that
%  normalize is only active if your custom protocol function calls
%  the function: [obsv] = norm_obsv( obsv, params, cfg ).
cfg.obsv_defs = { ...
  struct('name','XT', 'units','fraction',              'display',0, 'minplot',0, 'maxplot',1.05), ...
  struct('name','YP', 'units','fraction', 'norm','YT', 'display',0, 'minplot',0, 'maxplot',1.05), ...
  struct('name','YT', 'units','fraction',              'display',0, 'minplot',0, 'maxplot',1.05), ...
  struct('name','RP', 'units','fraction', 'norm','RT', 'display',1, 'minplot',0, 'maxplot',0.30), ...
  struct('name','RT', 'units','fraction',              'display',0, 'minplot',0, 'maxplot',1.05) ...
};
% initialize observable structs
cfg = init_observable_defs( cfg.obsv_defs, cfg );


% parameters for visualizing trajectories [REQUIRED only for visualization]
cfg.time_units = 's';                               % time units
cfg.sim_tstart = 0;                                 % simulation start time, s
cfg.sim_tstop  = 63;                                % simulation stop time, s
cfg.sim_dt     = 0.5;                               % time step for trajectory, s

% penalty for long integration times [OPTIONAL: useful for avoiding stiff parameter regions]
cfg.timepenalty = 0;


%% Parallel temperating options
% Defaults are usually ok. Things you may want to change: jobname, nchains,
% parallel, nswaps, adapt_last, energy_init_max, relstep_init.
% See core/init_config_defaults.m for a complete list of config options.
cfg.jobname = 'nfbkosc_pt';            % job name, for file input/output
cfg.parallel = 0;                      % parallel? true/false
cfg.maxlabs  = 4;                      % maximum number of labs for parallel computing
cfg.nchains  = 4;                      % number of chains
cfg.nswaps = 15000;                    % number of chain swaps (=number of saved samples!)
cfg.nsteps = 25;                       % number of steps between chain swaps
cfg.display_status_interval = 1;       % How often to display info
cfg.save_progress_interval = 1000;     % How often to save progress 
cfg.adapt_last = 4900;                 % last adaption step
cfg.energy_init_max = 100;             % maximum allowed energy for initialization
cfg.beta_init = 0.666;                 % beta initialization parameter
cfg.relstep_init = 0.05;               % relstep initialization parameter


%% Load experimental data [REQUIRED].
% The data file should import a cell array called "expt".
% Each cell array contains data for an experiment. It should define
%  the fields 'time', 'mean', 'stdev', 'nsamples' and 'weight'.
%  each array should have dimensions T x O, where T is the number of
%  time points and O is the number of observables, except for 'time'
%  which is is a column vector of length T. Set weight=0 if an observation
%  is missing or hidden (e.g. for future validation).  Missing data may be
%  indicated by NaN (not-a-number).
load data;   
% get number of experiments
cfg.nexpt = length(expt);
for d = 1 : cfg.nexpt
    % input signal (this is specific to the NFO model)
    cfg.data{d}.S  = expt{d}.S;
    % required fields for evaluating fit
    cfg.data{d}.time = expt{d}.time;
    cfg.data{d}.mean = expt{d}.mean;
    cfg.data{d}.stdev = expt{d}.stdev;
    cfg.data{d}.nsamples = expt{d}.nsamples;
    cfg.data{d}.weight   = expt{d}.weight;
end


%% load default functions
cfg = setup_default_functions( cfg );


%% setup custom function handles

% Sample parameter prior function [optional]:
%   prototype: [params] = @()
%
%cfg.sample_prior_fcn = @() ( <insert custom function here> );


% Log-prior pdf function [optional]:
%   prototype: [logpdf] = @(params)
%
%cfg.logpdf_prior_fcn = @(params) logpdf_prior( <insert custom function here> );


% Equilibration protocol function [optional]:
%   prototype: [err,sim,obsv] = @(params)
%
% NOTE: Equilibration is performed once per parameter set. If equilibration
% is dependent on the protocol, include equilibration in protocol function.
%
%cfg.equilibrate_fcn = @(params) ( <insert custom function here> );


% Define a simulation protocol for each experiment [REQUIRED]
%
% pass extra options to the protocol fcn here:
args = struct( ...
    'param_map',    cfg.param_map,    ...% useful for finding params by name
    'obsv_to_norm', cfg.obsv_to_norm, ...% required by "norm_obsv"
    'obsv_norm_by', cfg.obsv_norm_by  ...% required by "norm_obsv"
);
% experiment protocols here:
for d = 1 : cfg.nexpt
    % protocol specific parameters:
    S = cfg.data{d}.S;  % input signal to NFO
    % Experimental protocol function:
    %   prototype: [err,sim,obsv] = @(t,init,params)
    cfg.data{d}.protocol_fcn = @(t,init,params) simulate(t,init,params,S,args);
end


% Energy function [optional]
%   prototype: [energy] = @(params)
%
% NOTE: default is log-t-distribution (handles small samples better than square residual)
%
%cfg.energy_fcn = @(params) ( <insert custom function here> )




