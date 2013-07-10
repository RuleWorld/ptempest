function [cfg] = init_config_defaults()
%INIT_CONFIG_DEFAULTS set default configuration options
%
%   [cfg] = init_config_defaults()

    cfg = struct();

    % parallel temperating options
    cfg.jobname = 'pt';                    % job name, for file input/output
    cfg.shuffle  = 1;                      % shuffle random number streams (seed by clock)
    cfg.parallel = 0;                      % parallel? true/false
    cfg.maxlabs  = 4;                      % maximum number of labs for parallel computing
    cfg.nchains  = 4;                      % number of chains
    cfg.nswaps = 10000;                    % number of chain swaps
    cfg.nsteps = 25;                       % number of steps between chain swaps
    cfg.display_status_interval = 5;       % How often to display info
    cfg.save_progress_interval = 1000;     % How often to save progress 
    cfg.adapt_relstep_interval = 100;      % How often to adapt relative step size
    cfg.adapt_relstep_rate = 0.20;         % relative step size adaption rate
    cfg.optimal_step_acceptance = 0.24;    % optimal rate of step acceptance
    cfg.adapt_beta_interval = 250;         % how often to adapt temperature gradient
    cfg.adapt_beta_rate = 0.04;            % beta adaption rate
    cfg.optimal_swap_acceptance = 0.24;    % optimal rate of swap acceptance
    cfg.adapt_last = 2900;                 % last adaption step
    cfg.min_adaption_factor = 0.80;        % minimum adaption factor
    cfg.max_adaption_factor = 1.25;        % maximum adaption factor
    cfg.energy_init_max = 1000;            % maximum allowed energy for initialization
    cfg.max_init_steps = 500;              % maximum attempts at initialization
    cfg.max_beta = 1.0;                    % maximum chain beta (inverse of minimum chain temperature)
    cfg.beta_init = 0.666;                 % beta initialization parameter
    cfg.relstep_init = 0.01;               % relstep initialization parameter
    cfg.big_energy = 1e29;                 % big energy value: return if integrator fails
    cfg.timepenalty = 0;                   % penalty for long integration times

    % config trajectory visualization
    cfg.time_units = 's';                  % time units
    cfg.sim_tstart = 0;                    % simulation start time, s
    cfg.sim_tstop  = 20;                   % simulation stop time, s
    cfg.sim_dt     = 0.5;                  % time step for trajectory, s

end
