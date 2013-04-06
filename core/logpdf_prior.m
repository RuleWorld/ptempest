function [logpdf] = logpdf_prior( params, cfg )
%LOGPDF_PRIOR Calculate log-prior for parameters
%
%   [logpdf] = logpdf_prior( params, cfg )

% this for-loop version seems to be faster than the "functional" style below
logpdf = 0;
P = cfg.nparams;
for p = 1 : P
    logpdf = logpdf + cfg.param_logpdf{p}(params(p));
end

% this version is slower..
%logpdf = sum( arrayfun(@(p)(cfg.param_logpdf{p}(params(p))), 1:cfg.nparams) );
