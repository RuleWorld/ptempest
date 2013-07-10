function [] = plot_distr( params, cfg, figidx )
%PLOT_DISTR Plot default trajectory distribution over parameter ensemble
%
%  [] = plot_distr( params, cfg )
%  [] = plot_distr( params, cfg, figidx )
%
%  where 'params' is a (S x P) array of parameters and 'cfg' is the
%  configuration struct, and 'figidx' is an optional figure index.
%
%  S = number of samples, P = number of parameters


% permute params if the samples are lie along the 3rd index
if (size(params,1)==1 & size(params,3)>1 )
    params = permute(params,[3 2 1]);
end

% set upper and lower quantiles for distribution fill
upperquant    = 0.978;
midupperquant = 0.842;
midlowerquant = 0.158;
lowerquant    = 0.022;
% number of horizontal subplots
n_horiz_subplots = 4;
% legend location
legend_location = 'NorthEast';


% get index of observables to display
if isfield(cfg, 'obsv_display')
    obsv_disp = find(cfg.obsv_display);
else
    obsv_disp = [1 : cfg.nobsv];
end



% define xlabel
if isfield(cfg,'time_units')
    xlabel_string = sprintf( 'time (%s)', cfg.time_units );
else
    xlabel_string = 'time';
end
% set up time vector
sim_t = [cfg.sim_tstart : cfg.sim_dt : cfg.sim_tstop]';

% define ylabels
ylabel_strings = {};
for o = obsv_disp
    ylabel_strings{o} = sprintf( '%s (%s)', cfg.obsv_names{o}, cfg.obsv_units{o} );
end

% define experiment labels
expt_strings = {};
for d = 1:length(cfg.data)
    if isfield( cfg.data{d}, 'name' )
        expt_strings{d} = cfg.data{d}.name;
    else
        expt_strings{d} = sprintf( 'expt%d', d );
    end
end


% set up linestyles, etc (supports up to three experiments)
linestyle  = {'-b','-r','g','-k'};
fillstyle = {'b','r','g','k'};
errbar_style = {'o','s','v','x'};
errbar_color = {[0 0 0.7], [0.7 0 0], [0 0.7 0], [0.1 0.1 0.1]};


% dimensions
T = length(sim_t);     % time points
N = length(obsv_disp); % display observables
O = cfg.nobsv;         % observables
D = length(cfg.data);  % experiments
S = size(params,1);    % param samples


% initialize arrays for data storage
for d = 1:D
    sim_obsv{d} = zeros(T,O,S);
end


% loop over parameter samples
for s = 1:S

    % equillibrate, if required
    if isfield(cfg, 'equilibrate_fcn')
        [err,state,~] = cfg.equilibrate_fcn( params(s,:) );
        if (err)
            fprintf(1,'error equilibrating model\n');
            return;
        end
        init = state(end,:);
    else
        init = [];
    end

    % loop over experiments
    for d = 1:D

        % simulate experiment
        [err, ~, obsv] = cfg.data{d}.protocol_fcn( sim_t, init, params(s,:) );
        if (err) return; end
        % normalize and save results
        sim_obsv{d}(:,:,s) = norm_obsv( obsv, params(s,:), cfg );

    end

    % show progress
    fprintf(1,'.');
    if mod(s, 72)==0
        fprintf(1,'\n');
    end

end
fprintf(1,'\n');


%% calculate mean and stdev for simulated data
for d = 1:D
    % eliminate NaN's
    sim_obsv{d}( find( isnan(sim_obsv{d}) ) ) = 0;
    % compute quantiles
    sim_mean{d} = quantile(sim_obsv{d}, 0.5, 3);
    sim_max{d} = quantile(sim_obsv{d}, 1.0, 3);
    sim_upper{d} = quantile(sim_obsv{d}, upperquant, 3);
    sim_midupper{d} = quantile(sim_obsv{d}, midupperquant, 3);
    sim_midlower{d} = quantile(sim_obsv{d}, midlowerquant, 3);
    sim_lower{d} = quantile(sim_obsv{d}, lowerquant, 3);
    sim_min{d} = quantile(sim_obsv{d}, 0.0, 3);
end



% fetch experimental data
expt = cfg.data;
% set up errorbars
for d = 1:D
    errbar_time{d} = expt{d}.time;
    if isfield(cfg,'transform_data_for_plot_fcn')
        errbar_mid{d}  = cfg.transform_data_for_plot_fcn( expt{d}.mean );
        errbar_high{d} =  cfg.transform_data_for_plot_fcn( expt{d}.mean + expt{d}.stdev ) ...
                            - cfg.transform_data_for_plot_fcn( expt{d}.mean );
        errbar_low{d}  = cfg.transform_data_for_plot_fcn( expt{d}.mean ) ...
                            - cfg.transform_data_for_plot_fcn( expt{d}.mean - expt{d}.stdev );
    else
        errbar_mid{d}  = expt{d}.mean;
        errbar_high{d} = expt{d}.stdev;
        errbar_low{d}  = expt{d}.stdev;
    end
end


% plot experimental data and the distribution of model trajectories
if exist('figidx')  fh = figure(figidx);
else  fh = figure;  end
set( fh, 'Color', [1 1 1] );
for n = 1:N
    
    % get observable index
    o = obsv_disp(n);

    subplot( ceil(N/n_horiz_subplots), min(n_horiz_subplots,N), n);
    hold off;
    for d = 1:D

        %objh = fill( [sim_t; flipud(sim_t)], [sim_lower{d}(:,o); flipud(sim_upper{d}(:,o))], ...
        %             fillstyle{d}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'EdgeAlpha', 0.4);
        %set( get(get(objh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        %hold on;
        objh = fill( [sim_t; flipud(sim_t)], [sim_midlower{d}(:,o); flipud(sim_midupper{d}(:,o))], ...
                     fillstyle{d}, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'EdgeAlpha', 0.8);
        set( get(get(objh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        hold on;
        plot( sim_t, sim_mean{d}(:,o), linestyle{d}, 'linewidth', 1.1 );
       
    end

    for d = 1:D
        % plot errorbars (if any)
        if ( any(~isnan(expt{d}.mean(:,o))) )
            ebh = errorbar( errbar_time{d}, errbar_mid{d}(:,o), ...
                           errbar_low{d}(:,o), errbar_high{d}(:,o), ...
                           errbar_style{d}, 'color', errbar_color{d}, 'linewidth', 1.1 );
            % increase width of errorbars
            if ( isfield(cfg,'plot_xlim') )
                adjust_errorbars( ebh, cfg.plot_xlim );
            else
                adjust_errorbars( ebh, [cfg.sim_tstart cfg.sim_tstop] );
            end
        end
    end

    % add legend
    if (n==N)
        lh = legend( expt_strings );
        set(lh,'Location',legend_location);
    end

    % set axes limits
    if ( isfield(cfg,'plot_xlim') )
        xlim( cfg.plot_xlim );
    else
        xlim( [cfg.sim_tstart cfg.sim_tstop] );
    end
    if ( isfield(cfg,'plot_ylims') )
        ylim( cfg.plot_ylims{o} );
    end

    % label axes
    xlabel( xlabel_string, 'fontSize',12, 'interpreter','none');
    ylabel( ylabel_strings{o}, 'fontSize',12, 'interpreter','none');

    hold off;
end


% end plot_distr function
end



function [] = adjust_errorbars( ebh, xlims )
% increase width of errorbars  
    ebh = get(ebh,'children');
    bardata = get(ebh(2),'Xdata');
    left_idx = sort( [4:9:length(bardata), 7:9:length(bardata)] );
    right_idx = sort( [5:9:length(bardata), 8:9:length(bardata)] );
    % compute centers and width
    centerdata = (bardata(left_idx) + bardata(right_idx))/2;
    width = (xlims(2)-xlims(1))/40;
    % adjust errorbars
    bardata(left_idx)  = centerdata - width;
    bardata(right_idx) = centerdata + width;
    set(ebh(2),'Xdata',bardata);
end

