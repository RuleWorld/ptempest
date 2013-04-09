function [obsv] = norm_obsv( obsv, params, cfg )
%NORM_OBSV normalize observables
%
%   [obsv] = norm_obsv( obsv, params, cfg )

    % compute fractions
    obsv(:, cfg.obsv_to_norm) = obsv(:, cfg.obsv_to_norm) ./ obsv(:, cfg.obsv_norm_by);

return
