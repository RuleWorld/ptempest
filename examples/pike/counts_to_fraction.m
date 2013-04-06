function [frac] = counts_to_fraction( obsv, params, cfg )
% convert absolute count observables into frac over totals

    % compute fractions
    frac = obsv ./ obsv(:, cfg.obsv_norms);

return
