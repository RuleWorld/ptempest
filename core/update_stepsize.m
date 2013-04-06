function [stepsize] = update_stepsize( relative_step_size, cfg)
% update the stepsize array
stepsize = repmat(relative_step_size, [1 cfg.nparams]).*repmat(cfg.param_interval, [cfg.nchains 1]);   
