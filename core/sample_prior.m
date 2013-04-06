function [params] = sample_prior( cfg )
%SAMPLE_PRIOR Sample initial parameters from prior distribution
%
%   [params] = sample_prior( cfg )

P = cfg.nparams;
params = zeros(1,P);
for p = 1 : P
    params(p) = cfg.param_sample{p}();
end

