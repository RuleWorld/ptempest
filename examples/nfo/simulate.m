function [err, sp, obsv] = simulate( t, init, params, S, cfg )
%SIMULATE Generate trajectory for negative feedback oscillator (NFO)
%
%  [err, sp, obsv] = simulate( t, init, params, S, cfg )
%
%  generate negative feedback oscillator trajectory, where 't' is a vector
%  of time points, with parameters 'params', with input signal 'S',
%  and additional configuration arguments 'cfg'.
%
%  If 'init' = [], then default initial conditions are used. Note that
%  initial conditions are with respect to the first time point, not
%  necessarily t=0.  
%
%  Time units are seconds. The state trajectories and observables 
%  have concentration units


if isempty(init)
    init = [];
end

% set input signal
params( cfg.param_map.S ) = S;
% run simulation with default initial cond.
[err, ~, sp, obsv] = model( t, init, params, 1 );
if (err)
    sp = [];
    obsv = [];
    return;
end


return;

