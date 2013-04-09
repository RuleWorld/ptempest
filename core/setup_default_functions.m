function [cfg] = setup_default_functions( cfg )
%SETUP_DEFAULT_FUNCTIONS setup default ptempest functions
%
%   [cfg] = setup_default_functions( cfg )

    % suffix and regex for progress and init files
    cfg.progress_suffix = 'progress';
    cfg.progress_regex  = sprintf( '%s_%s*.mat', cfg.jobname, cfg.progress_suffix );
    cfg.init_suffix = 'init';
    cfg.init_regex  = sprintf( '%s_%s*.mat', cfg.jobname, cfg.init_suffix );

    % x-axis limits for plotting experiments
    cfg.plot_xlim = [cfg.sim_tstart, cfg.sim_tstop];

    % Sample parameters from prior.
    % prototype: [params] = @()
    cfg.sample_prior_fcn = @() sample_prior(cfg);

    % Compute log-prior probability of parameters.
    % prototype: [logpdf] = @(params)
    cfg.logpdf_prior_fcn = @(params) logpdf_prior(params,cfg);

    % proposal generator, gaussian default
    %   template: [params] = @(params,epsilon)
    cfg.proposal_fcn = @(params,epsilon) params + epsilon.*randn(1,cfg.nparams);

    % update stepsize
    %   template: [stepsize] = @(relstep)
    cfg.update_stepsize_fcn = @(relstep) repmat(relstep, [1 cfg.nparams]).*repmat(cfg.param_scale, [cfg.nchains 1]);

end
