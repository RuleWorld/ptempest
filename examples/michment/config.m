function [cfg] = config()
%CONFIG configure Michaelis-Menten model and parallel tempering options
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
cfg.param_defs = { ...
  struct('name','S0', 'prior','point', 'value', 1e2, 'units','nM'), ...
  struct('name','E0', 'prior','point', 'value', 1e1, 'units','nM'), ...
  struct('name','log_k1',  'prior','uniform', 'min', log10(1e-5), 'max',log10(1e1), 'units','log10 /nM/s'), ...
  struct('name','log_k2',  'prior','uniform', 'min', log10(1e-3), 'max',log10(1e3), 'units','log10 /s'  ), ...
  struct('name','log_k3',  'prior','uniform', 'min', log10(1e-3), 'max',log10(1e3), 'units','log10 /s'  )  ...
};
% initialize parameter distributions
cfg = init_parameter_defs( cfg.param_defs, cfg );


%% Observable definitions [REQUIRED]
%  Set display=0 to disable observable plot during visualization. Set
%  minplot/maxplot to see y-axis bounds in visualization scripts.
%  The 'norm' field should be the name of the observable that is used to
%  normalize the observable (e.g. YP is normalized by YT, the total quantity
%  of Y). Leave this field empty if normalization is not desired.
cfg.obsv_defs = { ...
  struct('name','S',  'units','nM', 'display',1, 'minplot',0, 'maxplot',1.1e2), ...
  struct('name','E',  'units','nM', 'display',1, 'minplot',0, 'maxplot',1.1e1), ...
  struct('name','ES', 'units','nM', 'display',1, 'minplot',0, 'maxplot',1.1e1), ...
  struct('name','P',  'units','nM', 'display',1, 'minplot',0, 'maxplot',1.1e2)
};
% initialize observable structs
cfg = init_observable_defs( cfg.obsv_defs, cfg );


% parameters for visualizing trajectories [REQUIRED only for visualization]
cfg.time_units = 's';                               % time units
cfg.sim_tstart = 0;                                 % simulation start time, s
cfg.sim_tstop  = 1010;                              % simulation stop time, s
cfg.sim_dt     = 5;                                 % time step for trajectory, s

% additional parameters required for simulation
cfg.S_species_idx = 1;        % index of S species (for setting concentrations)
cfg.E_species_idx = 2;        % index of E species (for setting concentrations)

% penalty for long integration times [OPTIONAL: useful for avoiding stiff parameter regions]
cfg.timepenalty = 0;


%% Parallel temperating options
% Defaults are usually ok. Things you may want to change: jobname, nchains,
% parallel, nswaps, adapt_last, energy_init_max, relstep_init.
% See core/init_config_defaults.m for a complete list of config options.
cfg.jobname = 'mm_pt';                 % job name, for file input/output
cfg.parallel = 0;                      % parallel? true/false
cfg.maxlabs  = 3;                      % maximum number of labs for parallel computing
cfg.nchains  = 3;                      % number of chains
cfg.nswaps = 8000;                     % number of chain swaps (=number of saved samples!)
cfg.nsteps = 25;                       % number of steps between chain swaps
cfg.display_status_interval = 10;      % How often to display info
cfg.save_progress_interval = 1000;     % How often to save progress 
cfg.adapt_last = 2900;                 % last adaption step
cfg.energy_init_max = 160;             % maximum allowed energy for initialization
cfg.beta_init = 0.5;                   % beta initialization parameter
cfg.relstep_init = 0.01;               % relstep initialization parameter


%% Load experimental data [REQUIRED].
% The data file should import a cell array called "expt".
% Each cell array contains data for an experiment. It should define
%  the fields 'time', 'mean', 'stdev', 'nsamples' and 'weight'.
%  each array should have dimensions T x O, where T is the number of
%  time points and O is the number of observables, except for 'time'
%  which is is a column vector of length T. Set weight=0 if an observation
%  is missing or hidden (e.g. for future validation).  Missing data may be
%  indicated by NaN (not-a-number).
load MM;   
% get number of experiments
cfg.nexpt = length(expt);
for d = 1 : cfg.nexpt
    % input signal (this is specific to the NFO model)
    cfg.data{d}.S0  = expt{d}.S0;
    cfg.data{d}.E0  = expt{d}.E0;
    % required fields for evaluating fit
    cfg.data{d}.time = expt{d}.time;
    cfg.data{d}.mean = expt{d}.mean;
    cfg.data{d}.stdev = expt{d}.stdev;
    %cfg.data{d}.nsamples = expt{d}.nsamples;
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


% Define a simulation protocol for each experiment [REQUIRED!]
%
% pass extra options to the protocol fcn here:
args = struct( ...
    'param_map',     cfg.param_map,     ...% useful for finding params by name
    'S_species_idx', cfg.S_species_idx, ...
    'E_species_idx', cfg.E_species_idx, ...
    'obsv_to_norm',  cfg.obsv_to_norm,  ...% required by "norm_obsv"
    'obsv_norm_by',  cfg.obsv_norm_by   ...% required by "norm_obsv"
);
% experiment protocols here:
for d = 1 : cfg.nexpt
    % protocol specific parameters:
    S0 = cfg.data{d}.S0;  % initial substrate
    E0 = cfg.data{d}.E0;  % initial substrate
    % Experimental protocol function:
    %   prototype: [err,sim,obsv] = @(t,init,params)
    cfg.data{d}.protocol_fcn = @(t,init,params) simulate(t,init,params,S0,E0,args);
end


% Energy function [REQUIRED!]
%   prototype: [energy] = @(params)
%
cfg.energy_fcn = @(params) energy_gaussian(params, cfg);




