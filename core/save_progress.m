function [] = save_progress( cfg, swap_idx, params_chain, energy_chain, step_acceptance, swap_acceptance, ...
                              swap_time, last_swap, relative_step_size, beta, relstep_history, beta_history );
% save parallel tempering progress
fprintf(1,'------------------------------------------------------\n');
fprintf(1,' SAVING PROGRESS . . . ');
savefile = sprintf( '%s_%s%06d.mat', cfg.jobname, cfg.progress_suffix, swap_idx );
save( savefile, 'params_chain', 'energy_chain', 'cfg', ...
                 'step_acceptance','swap_acceptance', ...
                 'swap_time', 'last_swap', ...
                 'relative_step_size', 'beta', ...
                 'relstep_history', 'beta_history', '-v7.3' );
% delete old progress files
% delete any progress files
progress_files = dir( cfg.progress_regex );
for file_idx = 1: length(progress_files)
    file = progress_files(file_idx);
    if ( ~strcmp(file.name,savefile) )
        delete( file.name );
    end
end
fprintf(1,'done\n');
