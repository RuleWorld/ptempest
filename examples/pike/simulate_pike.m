function [err, sp, obsv] = simulate_pike( t, init, params, egf0, egfr0, cfg )
%SIMULATE_PIKE Generate trajectories for pike model
%
%  [err, sp, obsv] = simulate_egfr_althyp( t, init, params, egf0, egf0, cfg )
%
%  generate Pike model trajectory where 't' is a vector of time points,
%  with parameters 'params', starting from initial conditions 'egf0' and 
%  'egfr0' (scalar),and model configuration 'cfg'.
%
%  If 'init' = [], then default initial conditions are used. Note that
%  initial conditions are with respect to the first time point, not
%  necessarily t=0.  
%
%  Time units are seconds. The state and observables trajectories have units
%  of absolute counts.


if isempty(init)

    % set EGF and EGFR initial conditions via parameters
    params( cfg.param_map.conc_Egf_0 )   = egf0;
    params( cfg.param_map.count_Egfr_0 ) = egfr0;

    % run simulation
    [err, ~, sp, obsv] = model_pike( t, [], params, 1 );
    if (err)
        sp = [];
        obsv = [];
        return;
    end

else

    % set EGF and EGFR initial conditions directly
    init( cfg.Egf_species_idx ) = egf0 * cfg.counts_per_mol / params(cfg.param_map.cell_dens);
    init( cfg.Egfr_species_idx ) = egfr0;

    % run simulation
    [err, ~, sp, obsv] = model_pike( t, init, params, 1 );
    if (err)
        sp = [];
        obsv = [];
        return;
    end

end

return;

