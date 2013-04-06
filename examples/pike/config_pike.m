function [cfg] = config_pike()
% CONFIG_PIKE configure model and parallel tempering options
%
%   [cfg] = config_pike()
%


% define parameters
scale = 1;
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


% observable definitions
cfg.obsv_defs = { ...
  struct('name','LigFree',  'units','fraction', 'norm','LigTotal', 'display',1, 'minplot',0, 'maxplot',1.05), ...
  struct('name','LigTotal', 'units','fraction', 'norm','LigTotal', 'display',0, 'minplot',0, 'maxplot',1.05), ...
  struct('name','RecBound', 'units','fraction', 'norm','RecTotal', 'display',1, 'minplot',0, 'maxplot',1.05), ...
  struct('name','RecDimer', 'units','fraction', 'norm','RecTotal', 'display',1, 'minplot',0, 'maxplot',1.05), ...
  struct('name','RecTotal', 'units','fraction', 'norm','RecTotal', 'display',0, 'minplot',0, 'maxplot',1.05) ...
};
% initialize observable structs
cfg = init_observable_defs( cfg.obsv_defs, cfg );



% simulation parameters
cfg.time_units = 's';            % time units
cfg.time = [0 1 10 100 1000]';   % time vector for equilibration
cfg.Egf_species_idx  = 1;        % index of Egf  species (for setting concentrations)
cfg.Egfr_species_idx = 2;        % index of Egfr species (for setting concentrations)

% unit conversions
cfg.counts_per_mol  = 6.0221e23;  % counts/mol
% parameters for controller errors, etc
cfg.big_energy = 1e29;    % large energy value
cfg.tolerance  = 1e-12;   % tolerance factor for comparing things to zero

% penalty for long integration times
cfg.timepenalty = 25;

% parameters for plotting simulations
cfg.sim_tstart = 0;                                 % simulation start time, sec
cfg.sim_tstop  = 120;                               % simulation stop time, sec
cfg.sim_dt     = 0.5;                               % time step for trajectory, sec
cfg.plot_xlim = [cfg.sim_tstart, cfg.sim_tstop];    % x-axis limits for plotting experiments


% load experimental data
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


% parallel temperating options
cfg.jobname = 'pike_pt';               % job name, for file input/output
cfg.shuffle  = 1;                      % shuffle random number streams (seed by clock)
cfg.parallel = 0;                      % parallel? true/false
cfg.maxlabs  = 2;                      % maximum number of labs for parallel computing
cfg.nchains  = 2;                      % number of chains
cfg.nswaps = 5000;                     % number of chain swaps
cfg.nsteps = 25;                       % number of steps between chain swaps
cfg.display_status_interval = 1;       % How often to display info
cfg.save_progress_interval = 250;      % How often to save progress 
cfg.adapt_relstep_interval = 100;      % How often to adapt relative step size
cfg.adapt_relstep_rate = 0.24;         % relative step size adaption rate
cfg.optimal_step_acceptance = 0.23;    % optimal rate of step acceptance
cfg.adapt_beta_interval = 250;         % how often to adapt temperature gradient
cfg.adapt_beta_rate = 0.05;            % beta adaption rate
cfg.optimal_swap_acceptance = 0.23;    % optimal rate of swap acceptance
cfg.adapt_last = 1200;                 % last adaption step
cfg.min_adaption_factor = 0.80;        % minimum adaption factor
cfg.max_adaption_factor = 1.25;        % maximum adaption factor
cfg.energy_init_max = 1200;            % maximum allowed energy for initialization
cfg.max_beta = 1.0;                    % maximum chain beta (inverse of minimum chain temperature)
cfg.beta_init = 0.56;                  % beta initialization parameter
cfg.relstep_init = 0.02;               % relstep initialization parameter

% suffix and regex for progress and init files
cfg.progress_suffix = 'progress';
cfg.progress_regex  = sprintf( '%s_%s*.mat', cfg.jobname, cfg.progress_suffix );
cfg.init_suffix = 'init';
cfg.init_regex  = sprintf( '%s_%s*.mat', cfg.jobname, cfg.init_suffix );


%% define custom function handles
% initial proposal generator
%   template: [params] = @()
cfg.sample_prior_fcn = @() sample_prior(cfg);

% parameter prior
%   template: [logpdf] = @(params)
cfg.logpdf_prior_fcn = @(params) logpdf_prior(params,cfg);

% equilibration protocol
%   template: [err,sim,obsv] = @(params)
%cfg.equilibrate_fcn = @(params) ( ... );

% simulation protocols
%   template: [err,sim,obsv] = @(t,init,params)
% gather arguments required by protocols
args = struct( ...
    'param_map', cfg.param_map, ...
    'counts_per_mol', cfg.counts_per_mol, ...
    'Egf_species_idx',  cfg.Egf_species_idx, ...
    'Egfr_species_idx', cfg.Egfr_species_idx, ...
    'obsv_norms', cfg.obsv_norms ...
);
for d = 1 : cfg.nexpt
    egf0  = cfg.data{d}.egf;
    egfr0 = cfg.data{d}.egfr;
    cfg.data{d}.protocol_fcn = @(t,init,params) simulate_pike(t,init,params,egf0,egfr0,args);
end

% energy function, use "tdistr" if number of samples is small
%   template: [energy] = @(params)
cfg.energy_fcn = @(params) energy_tdistr(params,cfg);

% proposal generator, gaussian default
%   template: [params] = @(params,epsilon)
cfg.proposal_fcn = @(params,epsilon) params + epsilon.*randn(1,cfg.nparams);

% update stepsize
%   template: [stepsize] = @(relstep)
cfg.update_stepsize_fcn = @(relstep) repmat(relstep, [1 cfg.nparams]).*repmat(cfg.param_scale, [cfg.nchains 1]);



