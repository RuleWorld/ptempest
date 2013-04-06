function [cfg] = init_observable_defs( obsv_defs, cfg )
%INIT_OBSERVABLE_DEFS read observable definitions and set up structs

    cfg.nobsv = length(obsv_defs);
    cfg.obsv_names = {};
    cfg.obsv_units = {};
    cfg.obsv_display = ones(1,cfg.nobsv);
    cfg.plot_ylims = {};

    found_norm  = 0;
    for o = 1:cfg.nobsv
        % map from observable names to indices
        cfg.obsv_map.(obsv_defs{o}.name) = o;
        % set up structs
        cfg.obsv_names{o} = obsv_defs{o}.name;
        cfg.obsv_units{o} = obsv_defs{o}.units;
        if isfield(obsv_defs{o}, 'norm')
            found_norm = 1;            
        end
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

    if found_norm
        % find index of observable for normalization
        cfg.obsv_norms = zeros(1,cfg.nobsv);
        for o = 1:cfg.nobsv
            cfg.obsv_norms(o) = cfg.obsv_map.(obsv_defs{o}.norm);
        end
    end

end
