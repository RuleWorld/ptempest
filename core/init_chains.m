function [energy_curr, params_curr] = init_chains( epsilon, cfg )
% find initial starting point
fprintf(1,'------------------------------------------------------\n');
fprintf(1,'Looking for good initial parameter set . . .\n');
energy_init = Inf;
params_init = zeros(1, cfg.nparams);
counter = 1;
for counter = 1 : cfg.max_init_steps
    params = cfg.sample_prior_fcn();
    energy = cfg.energy_fcn( params );
    if (energy < energy_init)
        energy_init = energy;
        params_init = params;
    end
    if ( mod(counter,5)==0 )
        fprintf( 1, '  . . . %d attempts (energy curr=%g, best=%g)\n', counter, energy, energy_init );
    end
    if ( energy_init < cfg.energy_init_max )
        break;
    end
end
fprintf(1,' . . . energy at initial starting point is %g\n', energy_init );
% initialize chains
fprintf(1,'------------------------------------------------------\n');
fprintf(1,'Initializing chains...\n');
params_curr = zeros(cfg.nchains,cfg.nparams);
energy_curr = zeros(cfg.nchains,1);
for chain_idx = 1 : cfg.nchains
    fprintf(1,'chain %d . . .  ', chain_idx );
    energy_curr(chain_idx) = Inf;
    for counter = 1 : cfg.max_init_steps
        params = cfg.proposal_fcn( params_init, epsilon(chain_idx,:) ); 
        energy = cfg.energy_fcn( params );
        if ( energy < energy_curr(chain_idx) )
            params_curr(chain_idx,:) = params; 
            energy_curr(chain_idx) = energy;
        end
        if ( energy_curr(chain_idx) < cfg.energy_init_max )
            break;
        end            
    end
    fprintf(1,'initial energy = %g\n', energy_curr(chain_idx) );
end
