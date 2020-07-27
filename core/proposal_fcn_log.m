function [pout] = proposal_fcn_log(params, epsilon, cfg)
% proposal_fcn_log -- Modified proposal function that handles log scaled
% parameters.
% It first transforms loguniform parameters specified by logparams argument
% into log space, then translates all parameters, then transforms log params
% back to original space.

if (~exist('cfg.logparams'))
    lp=boolean([]);
    for i=1:length(cfg.param_defs)
        lp(end+1)= strcmp(cfg.param_defs{i}.prior,'loguniform');
    end
    cfg.logparams=lp;
end
logparams=cfg.logparams;
pout=params;
% Transform logscale parameters
pout(logparams)= log10(pout(logparams));

% Translate parameters using Gaussian proposal
pout= pout + epsilon.*randn(1,cfg.nparams);

% Back transform logscale parameters
pout(logparams)= 10.^(pout(logparams));

end

