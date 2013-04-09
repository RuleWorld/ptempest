function [energy] = energy_tdistr( params, cfg )
%ENERGY_TDISTR Calculate parameter energy using t-distribution likelihood model
% 
%  [energy] = energy_tdistr( params, cfg )

% check prior probability on parameters
logprior = cfg.logpdf_prior_fcn(params);
if isinf(logprior)
    energy = cfg.big_energy;
    return;
end
energy = -logprior;

% start integration timer
simtimer = tic;

% equillibrate, if required
if isfield(cfg, 'equilibrate_fcn')
    [err,state,~] = cfg.equilibrate_fcn( params );
    if (err)
        energy = cfg.big_energy;
        return;
    end
    init = state(end,:);
else
    init = [];
end

% simulate experiments
for d = 1:cfg.nexpt

    % simulate experiment (default initial conditions)
    [err,~,obsv] = cfg.data{d}.protocol_fcn( cfg.data{d}.time, init, params );
    if (err)
        energy = cfg.big_energy;
        return;
    end

    % normalize obsv
    [obsv] = norm_obsv( obsv, params, cfg );

    % heuristic penalties (do this before transformating obsv)
    if isfield(cfg.data{d}, 'heuristic_penalty_fcn')
        penalty = cfg.data{d}.heuristic_penalty_fcn(obsv, params);
        if isinf(penalty)
            energy = cfg.big_energy;
            return;
        end
        energy = energy + penalty;
    end

    % if necessary, transform simulated trajectory for computing fitness
    if isfield(cfg, 'transform_sim_for_fit_fcn')
        obsv = cfg.transform_sim_for_fit_fcn(obsv,params);
    end

    % compute t-statistics
    t = (cfg.data{d}.mean - obsv)./(cfg.data{d}.stdev./sqrt(cfg.data{d}.nsamples));
    loglike = nansum(nansum( cfg.data{d}.weight .* tlogpdf(t, cfg.data{d}.nsamples) ));

    % subtract likelihood from energy
    energy = energy - loglike;

end

% penalize for slow integrations
dt = toc(simtimer);
energy = energy + cfg.timepenalty*dt^2;

% all done
return;
