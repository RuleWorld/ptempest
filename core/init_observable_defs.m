function [cfg] = init_observable_defs( obsv_defs, cfg )
%INIT_OBSERVABLE_DEFS read observable definitions and set up structs

    cfg.nobsv = length(obsv_defs);
    cfg.obsv_names = {};
    cfg.obsv_units = {};
    cfg.obsv_display = ones(1,cfg.nobsv);
    cfg.plot_ylims = {};
    cfg.obsv_to_norm = [];  % index of observables that will be normalized
    cfg.obsv_norm_by = [];  % index of observables used to normalize the above

    for o = 1:cfg.nobsv
        % map from observable names to indices
        cfg.obsv_map.(obsv_defs{o}.name) = o;
        % set up structs
        cfg.obsv_names{o} = obsv_defs{o}.name;
        cfg.obsv_units{o} = obsv_defs{o}.units;
        if isfield(obsv_defs{o}, 'display')
            cfg.obsv_display(o) = obsv_defs{o}.display;
        end
        % parse max/minplot fields
        if isfield(obsv_defs{o}, 'maxplot') 
            if isfield(obsv_defs{o}, 'minplot')
                cfg.plot_ylims{o} = [obsv_defs{o}.minplot obsv_defs{o}.maxplot];
            else
                cfg.plot_ylims{o} = [0 obsv_defs{o}.maxplot];
            end
        elseif isfield(obsv_defs{o}, 'minplot') 
            cfg.plot_ylims{o} = [obsv_defs{o}.minplot Inf];
        else
            cfg.plot_ylims{o} = [0 Inf];
        end
    end

    % do a second pass to find normalizations
    for o = 1:cfg.nobsv
        if isfield(obsv_defs{o}, 'norm')
            cfg.obsv_to_norm = [cfg.obsv_to_norm, o];
            cfg.obsv_norm_by = [cfg.obsv_norm_by, cfg.obsv_map.(obsv_defs{o}.norm) ];
        end
    end

end
