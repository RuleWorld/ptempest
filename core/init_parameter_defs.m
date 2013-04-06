function [cfg] = init_parameter_defs( param_defs, cfg )
%INIT_PARAMETER_DEFS read parameter definitions and set up distributions
    cfg.nparams = length(cfg.param_defs);
    cfg.param_names = {};
    cfg.param_sample = {};
    cfg.param_logpdf = {};
    cfg.param_scale = zeros(1,cfg.nparams);
    cfg.param_location = zeros(1,cfg.nparams);
    for p = 1:cfg.nparams
        param_def = cfg.param_defs{p};
        % map from parameter names to indices
        cfg.param_map.(param_def.name) = p;
        % set up distribution
        cfg.param_names{p} = param_def.name;
        if (strcmp(param_def.prior,'point'))
            cfg.param_sample{p} = @() pointrnd( param_def.value );
            cfg.param_logpdf{p} = @(x) pointlogpdf( x, param_def.value );
            cfg.param_pdf{p} = @(x) pointpdf( x, param_def.value );
            cfg.param_location(p) = param_def.value;
            cfg.param_scale(p) = 0;
        elseif (strcmp(param_def.prior,'uniform'))
            cfg.param_sample{p} = @() unifrnd( param_def.min, param_def.max );
            cfg.param_logpdf{p} = @(x) uniflogpdf( x, param_def.min, param_def.max );
            cfg.param_pdf{p} = @(x) unifpdf( x, param_def.min, param_def.max );
            cfg.param_location(p) = (param_def.min + param_def.max)/2;
            cfg.param_scale(p) = (param_def.max - param_def.min)/2;
        elseif (strcmp(param_def.prior,'normal'))
            cfg.param_sample{p} = @() normrnd( param_def.mu, param_def.sigma );
            cfg.param_logpdf{p} = @(x) normlogpdf( x, param_def.mu, param_def.sigma );
            cfg.param_pdf{p} = @(x) normpdf( x, param_def.mu, param_def.sigma );
            cfg.param_location(p) = param_def.mu;
            cfg.param_scale(p) = param_def.sigma;
        elseif (strcmp(param_def.prior,'lognormal'))
            cfg.param_sample{p} = @() lognrnd( param_def.mu, param_def.sigma );
            cfg.param_logpdf{p} = @(x) lognlogpdf( x, param_def.mu, param_def.sigma );
            cfg.param_pdf{p} = @(x) lognpdf( x, param_def.mu, param_def.sigma );
            cfg.param_location(p) = exp(param_def.mu + param_def.sigma^2/2);
            cfg.param_scale(p) = sqrt((exp(param_def.sigma^2) - 1).*exp(2*param_def.mu + param_def.sigma^2));
        elseif (strcmp(param_def.prior,'laplace'))
            cfg.param_sample{p} = @() laplacernd( param_def.mu, param_def.b );
            cfg.param_logpdf{p} = @(x) laplacelogpdf( x, param_def.mu, param_def.b );
            cfg.param_pdf{p} = @(x) laplacepdf( x, param_def.mu, param_def.b );
            cfg.param_location(p) = param_def.mu;
            cfg.param_scale(p) = param_def.b;
        elseif (strcmp(param_def.prior,'boundedlaplace'))
            cfg.param_sample{p} = @() boundedlaplacernd( param_def.mu, param_def.b, param_def.min, param_def.max );
            cfg.param_logpdf{p} = @(x) boundedlaplacelogpdf( x, param_def.mu, param_def.b, param_def.min, param_def.max );
            cfg.param_pdf{p} = @(x) boundedlaplacepdf( x, param_def.mu, param_def.b,param_def.min, param_def.max );
            cfg.param_location(p) = param_def.mu;
            cfg.param_scale(p) = param_def.b;
        elseif (strcmp(param_def.prior,'beta'))
            cfg.param_sample{p} = @() betarnd( param_def.alpha, param_def.beta );
            cfg.param_logpdf{p} = @(x) betalogpdf( x, param_def.alpha, param_def.beta );
            cfg.param_pdf{p} = @(x) betapdf( x, param_def.alpha, param_def.beta );
            cfg.param_location(p) = 0.5;
            cfg.param_scale(p) = 0.5;
        elseif (strcmp(param_def.prior,'exponential'))
            cfg.param_sample{p} = @() exprnd( param_def.mu );
            cfg.param_logpdf{p} = @(x) explogpdf( x, param_def.mu );
            cfg.param_pdf{p} = @(x) exppdf( x, param_def.mu );
            cfg.param_location(p) = param_def.mu;
            cfg.param_scale(p) = param_def.mu;
        else
            fprintf(1, 'Parameter prior distribution %s is not recognized (parameter %s), try again.\n', ...
                param_def.prior, param_def.name );
        end
    end

end
