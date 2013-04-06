function [energy_curr, params_curr] = init_chains( epsilon, cfg )
% find initial starting point
fprintf(1,'------------------------------------------------------\n');
fprintf(1,'Looking for good initial parameter set . . .\n');
energy_init = Inf;
min_energy_sofar = Inf;
params_init = zeros(1, cfg.nparams);
counter = 1;
while ( energy_init > cfg.energy_init_max )
    params_init = cfg.sample_prior_fcn();
    energy_init = cfg.energy_fcn( params_init );
    min_energy_sofar = min( [energy_init, min_energy_sofar] );
    counter = counter + 1;
    if ( mod(counter,1)==0 )
        fprintf( 1, '  . . . %d attempts (energy curr=%g, best=%g)\n', counter, energy_init, min_energy_sofar );
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
    while ( energy_curr(chain_idx) > cfg.energy_init_max )
        params_curr(chain_idx,:) = cfg.proposal_fcn( params_init, epsilon(chain_idx,:) ); 
        energy_curr(chain_idx) = cfg.energy_fcn( params_curr(chain_idx,:) );
    end
    fprintf(1,'initial energy = %g\n', energy_curr(chain_idx) );
end
