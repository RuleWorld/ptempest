function [cfg] = config_pike()
% CONFIG_PIKE configure model and parallel tempering options
%
%   [cfg] = config_pike()
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
scale = 1;       % global weight parameter for prior strength (or "width"). scale>0.
RT = 2.5774863;
% parameter definitions (see BNGL model for parameter descriptions)
cfg.param_defs = { ...
  struct('name','RT', 'prior','point', 'value', 2.5774863, 'units','kJ/mol'), ...
  struct('name','NA', 'prior','point', 'value', 6.0221e23, 'units','/mol'  ), ...
  struct('name','PI', 'prior','point', 'value', 3.1415927, 'units','(none)'), ...
  struct('name','radcell',   'prior','point', 'value', 1e-4, 'units','dm'    ), ...
  struct('name','cell_dens', 'prior','point', 'value', 1e8,  'units','/L'    ), ...
  struct('name','width_PM',  'prior','point', 'value', 1e-7, 'units','dm'    ), ...
  struct('name','conc_Egf_0',   'prior','point', 'value',20e-9, 'units','M'),     ...
  struct('name','count_Egfr_0', 'prior','point', 'value',24000, 'units','/cell'), ...
  struct('name','Gf_LR',   'prior','laplace',  'mu',-57.35, 'b',scale*RT*log(2), 'units','kJ/mol'), ...
  struct('name','Gf_RR',   'prior','laplace',  'mu',-29.82, 'b',scale*RT*log(2), 'units','kJ/mol'), ...
  struct('name','Gf_LRR',  'prior','laplace',  'mu',-0.365, 'b',scale*RT*log(2), 'units','kJ/mol'), ...
  struct('name','Gf_LRRL', 'prior','laplace',  'mu',7.079,  'b',scale*RT*log(2), 'units','kJ/mol'), ...
  struct('name','phi',   'prior','point', 'value',0.5,   'units','(none)'), ...
  struct('name','Ea_LR', 'prior','point', 'value',-17.8, 'units','kJ/mol'), ...
  struct('name','Ea_RR', 'prior','point', 'value',-11.9, 'units','kJ/mol')  ...
};
% initialize parameter distributions
cfg = init_parameter_defs( cfg.param_defs, cfg );


%% Observable definitions [REQUIRED]
%  Set display=0 to disable observable plot during visualization. Set
%  minplot/maxplot to see y-axis bounds in visualization scripts.
%  The 'norm' field should be the name of the observable that is used to
%  normalize the observable. Leave this field empty if normalization is not
%  desired.
cfg.obsv_defs = { ...
  struct('name','LigFree',  'units','fraction', 'norm','LigTotal', 'display',1, 'minplot',0, 'maxplot',1.05), ...
  struct('name','LigTotal', 'units','fraction',                    'display',0, 'minplot',0, 'maxplot',1.05), ...
  struct('name','RecBound', 'units','fraction', 'norm','RecTotal', 'display',1, 'minplot',0, 'maxplot',1.05), ...
  struct('name','RecDimer', 'units','fraction', 'norm','RecTotal', 'display',1, 'minplot',0, 'maxplot',1.05), ...
  struct('name','RecTotal', 'units','fraction',                    'display',0, 'minplot',0, 'maxplot',1.05) ...
};
% initialize observable structs
cfg = init_observable_defs( cfg.obsv_defs, cfg );


% parameters for visualizing trajectories [REQUIRED only for visualization]
cfg.time_units = 's';                               % time units
cfg.sim_tstart = 0;                                 % simulation start time, sec
cfg.sim_tstop  = 120;                               % simulation stop time, sec
cfg.sim_dt     = 0.5;                               % time step for trajectory, sec

% penalty for long integration times
cfg.timepenalty = 25;

% additional simulation parameters
cfg.time = [0 1 10 100 1000]';   % time vector for equilibration
cfg.Egf_species_idx  = 1;        % index of Egf  species (for setting concentrations)
cfg.Egfr_species_idx = 2;        % index of Egfr species (for setting concentrations)

% unit conversions
cfg.counts_per_mol  = 6.0221e23;  % counts/mol


%% Parallel temperating options
% Defaults are usually ok. Things you may want to change: jobname, nchains,
% parallel, nswaps, adapt_last, energy_init_max, relstep_init.
% See core/init_config_defaults.m for a complete list of config options.
cfg.jobname = 'pike_pt';               % job name, for file input/output
cfg.parallel = 0;                      % parallel? true/false
cfg.maxlabs  = 3;                      % maximum number of labs for parallel computing
cfg.nchains  = 3;                      % number of chains
cfg.nswaps = 8000;                     % number of chain swaps (=number of saved samples!)
cfg.nsteps = 25;                       % number of steps between chain swaps
cfg.display_status_interval = 5;       % How often to display info
cfg.save_progress_interval = 1000;     % How often to save progress 
cfg.adapt_last = 2900;                 % last adaption step
cfg.energy_init_max = 180;             % maximum allowed energy for initialization
cfg.beta_init = 0.56;                  % beta initialization parameter
cfg.relstep_init = 0.02;               % relstep initialization parameter


%% Load experimental data [REQUIRED].
% The data file should import a cell array called "expt".
% Each cell array contains data for an experiment. It should define
%  the fields 'time', 'mean', 'stdev', 'nsamples' and 'weight'.
%  each array should have dimensions T x O, where T is the number of
%  time points and O is the number of observables, except for 'time'
%  which is is a column vector of length T. Set weight=0 if an observation
%  is missing or hidden (e.g. for future validation).  Missing data may be
%  indicated by NaN (not-a-number).
load synthdata_pike;
% get number of experiments
cfg.nexpt = length(expt);
for d = 1 : cfg.nexpt
    % EGF dose and number of EGFR receptors
    cfg.data{d}.egf  = expt{d}.egf;
    cfg.data{d}.egfr = expt{d}.egfr;
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


% Define a simulation protocol for each experiment [REQUIRED!]
%
% pass extra options to the protocol fcn here:
args = struct( ...
    'param_map',        cfg.param_map,        ...% useful for finding params by name
    'counts_per_mol',   cfg.counts_per_mol,   ...
    'Egf_species_idx',  cfg.Egf_species_idx,  ...
    'Egfr_species_idx', cfg.Egfr_species_idx, ...
    'obsv_to_norm',     cfg.obsv_to_norm,     ...% required by "norm_obsv"
    'obsv_norm_by',     cfg.obsv_norm_by      ...% required by "norm_obsv"
);
% experiment protocols here:
for d = 1 : cfg.nexpt
    % protocol specific parameters:
    egf0  = cfg.data{d}.egf;
    egfr0 = cfg.data{d}.egfr;
    % Experimental protocol function:
    %   prototype: [err,sim,obsv] = @(t,init,params)
    cfg.data{d}.protocol_fcn = @(t,init,params) simulate_pike(t,init,params,egf0,egfr0,args);
end


% Energy function [REQUIRED!]
%   prototype: [energy] = @(params)
%
cfg.energy_fcn = @(params) energy_tdistr(params, cfg);


