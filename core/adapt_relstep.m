function [relative_step_size] = adapt_relstep( relative_step_size, step_acceptance, swap_idx, cfg )
% adapt relative step size
fprintf(1,'------------------------------------------------------\n');
fprintf(1,' Adapting relative step size . . .\n');
% based on acceptance rate of all chains
step_acceptance_rate = sum( step_acceptance(:, (swap_idx-cfg.adapt_relstep_interval+1):swap_idx), 2 ) ...
                        / (cfg.nsteps*cfg.adapt_relstep_interval);
for chain_idx = 1 : cfg.nchains                                                                 
    old_relstep = relative_step_size(chain_idx);
    adaption_factor = (step_acceptance_rate(chain_idx)/cfg.optimal_step_acceptance)^cfg.adapt_relstep_rate;
    adaption_factor = min( [max([cfg.min_adaption_factor, adaption_factor]), cfg.max_adaption_factor] );
    relative_step_size(chain_idx) = adaption_factor*old_relstep;
    fprintf(1,'  chain %d: recent acceptance=%-8.3g old relstep=%-8.3g new relstep=%-8.3g\n', ...
       chain_idx, step_acceptance_rate(chain_idx), old_relstep, relative_step_size(chain_idx) );
end    

