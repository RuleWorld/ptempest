function [] = save_init( cfg, params_chain, energy_chain, beta_history, relstep_history );
% save an 'init' file for continuing this chain
fprintf(1,'------------------------------------------------------\n');
fprintf(1,' CREATING INIT FILE (for chain continuation) . . . ');
% create init structure
init.params = params_chain(:,:,end);
init.energy = energy_chain(:,end);
init.beta   = beta_history(:,end);
init.relstep = relstep_history(:,end);
% save structure
savefile = sprintf( '%s_%s.mat', cfg.jobname, cfg.init_suffix );
save( savefile, 'init', '-v7.3' );
fprintf(1,'done\n');
