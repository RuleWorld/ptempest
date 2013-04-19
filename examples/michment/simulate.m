function [err, sp, obsv] = simulate( t, init, params, S0, E0, cfg )
%SIMULATE Generate trajectory for Michaelis-Menten
%
%  [err, sp, obsv] = simulate( t, init, params, S0, E0, cfg )
%
%  generate Michaelis-Menten trajectory, where 't' is a vector
%  of time points, with parameters 'params', with initial substrate 'S0',
%  initial enzyme 'S0', and additional configuration arguments 'cfg'.
%
%  If 'init' = [], then default initial conditions are used. Note that
%  initial conditions are with respect to the first time point, not
%  necessarily t=0.  
%
%  Time units are seconds. The state trajectories and observables 
%  have concentration units


if isempty(init)
    % set initial S and E using parameters
    init = [];
    params( cfg.param_map.S0 ) = S0;
    params( cfg.param_map.E0 ) = E0;

else
    % set S and E as initial conditions
    init( cfg.S_species_idx ) = S0;
    init( cfg.E_species_idx ) = E0;

end

% run simulation
[err, ~, sp, obsv] = model( t, init, params, 1 );
if (err)
    sp = [];
    obsv = [];
    return;
end

return;

