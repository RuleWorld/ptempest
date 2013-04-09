function [err] = plot_pike( params, cfg, figidx )
%PLOT_PIKE Simulate and plot parameter scan over ligand and receptor concentration
% 
%  [] = plot_pike( params, cfg, figidx )


% permute params if the samples are lie along the 3rd index
if (size(params,1)==1 & size(params,3)>1 )
    params = permute(params,[3 2 1]);
end

% set upper and lower quantiles for distribution fill
upperquant    = 0.978;
midupperquant = 0.842;
medquant = 0.5;
midlowerquant = 0.158;
lowerquant    = 0.022;
% number of horizontal subplots
n_horiz_subplots = 3;
% legend location
legend_location = 'NorthWest';

% EGF doses (M)
x_obsv_idx = 1;
egf = [          1e-11, 1.72102769e-11,  2.9619363e-11, 5.09757437e-11, 8.77306662e-11, ...
        1.50986905e-10, 2.59852645e-10, 4.47213595e-10, 7.69666979e-10, 1.32461818e-09, ...
        2.27970456e-09, 3.92343467e-09, 6.75233969e-09, 1.16209635e-08,          2e-08  ];
nL = length(egf);

% EGFR expression levels (receptors/cell)
egfr = [ 24000, 43000, 92000, 120000, 231000, 447000 ];
nR = length(egfr);


% get index of observables to display
%if isfield(cfg, 'obsv_display')
%    obsv_disp = find(cfg.obsv_display);
%else
%    obsv_disp = [1 : cfg.nobsv];
%end
obsv_disp = [5 7];


% dimensions
N = length(obsv_disp); % display observables
O = cfg.nobsv;         % observables
S = size(params,1);    % param samples

% x/y axis limits
if ( isfield(cfg,'plot_xlim') )
    %xlim( cfg.plot_xlim );
    plot_xlim = [7e-12 2.5e-8];
else
    plot_xlim = [7e-12 2.5e-8];
end

for o = 1:O
    if ( isfield(cfg,'plot_ylims') )
        %plot_ylims{o} = cfg.plot_ylims{o};
        plot_ylims{o} = [0 1.05];
    else
        plot_ylims{o} = [0 1.05];
    end
end


% style for errorbars
errbar_style = {'o','o','o','o','o','o'};
errbar_color = {[0 0 0.7], [0 0.7 0], [0.7 0 0], [0 0.7 0.7], [0.7 0 0.7],[0.2 0.2 0.2]};
% style for lines
line_style = {'-b','-g','-r','-c','-m','-k'};
% style for fill
fill_style = {'b','g','r','c','m','k'};

% xlabel
xlabel_string = 'log [Free EGF] (M)';

% define ylabels
ylabel_strings = {};
for o = obsv_disp
    ylabel_strings{o} = sprintf( '%s (%s)', cfg.obsv_names{o}, cfg.obsv_units{o} );
end

% experiment stings
expt_strings = { 'R/cell=24k', 'R/cell=43k', 'R/cell=92k', 'R/cell=120k', 'R/cell=231k', 'R/cell=447k' };

% simulate experiments
expt = {};
for r = 1 : nR

    expt{r}.egf  = egf;
    expt{r}.egfr = egfr(r);
    expt{r}.obsv = zeros(1,cfg.nobsv,S);
    expt{r}.obsv_lower = zeros(length(egf),cfg.nobsv);
    expt{r}.obsv_med   = zeros(length(egf),cfg.nobsv);
    expt{r}.obsv_upper = zeros(length(egf),cfg.nobsv);
    expt{r}.mean   = zeros(length(egf),cfg.nobsv);
    expt{r}.stderr = zeros(length(egf),cfg.nobsv);

    for l = 1 : nL

        % get experiment index
        n = (r-1)*nL + l;

        % loop over parameter sets
        for s = 1:S

            % set EGF and EGFR initial conditions via parameters
            par = params(s,:);
            par( cfg.param_map.conc_Egf_0 )   = egf(l);
            par( cfg.param_map.count_Egfr_0 ) = egfr(r);

            % run simulation
            [err, ~, sp, obsv] = model_pike( cfg.time, [], par, 1 );
            if (err)
                sp = [];
                obsv = [];
                return;
            end
        
            % transform observables from counts to fractions
            [obsv] = norm_obsv(obsv(end,:),params,cfg);

            % convert free ligand to concentration
            totals = [egf(l) 1 1 1 1 1 1 1 ];
            obsv = obsv .* totals;

            expt{r}.obsv(1,:,s) = obsv;
    
        end

        expt{r}.obsv_lower(l,:) = quantile( expt{r}.obsv(1,:,:), lowerquant, 3);
        expt{r}.obsv_med(l,:)   = quantile( expt{r}.obsv(1,:,:), medquant, 3);
        expt{r}.obsv_upper(l,:) = quantile( expt{r}.obsv(1,:,:), upperquant, 3);


        % get lab measurements
        expt{r}.mean(l,:)   = cfg.data{n}.mean(end,:) .* totals;
        expt{r}.stderr(l,:) = cfg.data{n}.stdev(end,:) .* totals ./ sqrt(cfg.data{n}.nsamples(end,:)) ;

    end

end

% open figure
if exist('figidx') fh = figure(figidx);
else fh = figure; end
set( fh, 'Color', [1 1 1] );
set( fh,'PaperPositionMode','auto');
set( fh,'PaperOrientation','landscape');
hold off;


for n = 1:N
    
    % get observable index
    o = obsv_disp(n);

    sph = subplot( ceil(N/n_horiz_subplots), min(n_horiz_subplots,N), n);
    for r = 1:nR

        objh = fill( log10([expt{r}.obsv_med(:,x_obsv_idx); flipud(expt{r}.obsv_med(:,x_obsv_idx))]), ...
                     [expt{r}.obsv_lower(:,o); flipud(expt{r}.obsv_upper(:,o)) ], ...
                     fill_style{r}, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.7 ); %'EdgeColor', 'none', 'EdgeAlpha', 0.8 
        set( get(get(objh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        hold on;
        %objh = fill( [sim_t; flipud(sim_t)], [sim_midlower{d}(:,o); flipud(sim_midupper{d}(:,o))], ...
        %             fillstyle{d}, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'EdgeAlpha', 0.8);
        %set( get(get(objh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        %hold on;
        plot( log10(expt{r}.obsv_med(:,x_obsv_idx)), expt{r}.obsv_med(:,o), line_style{r}, 'linewidth', 1.1 );
       
    end

    for r = 1:nR
        % plot errorbars (if any)
        if (n~=1) continue; end;
        if ( any(~isnan(expt{r}.mean(:,o))) )
            ebh = errorbar( log10(expt{r}.mean(:,x_obsv_idx)), expt{r}.mean(:,o), ...
                            expt{r}.stderr(:,o), expt{r}.stderr(:,o), ...
                            errbar_style{r}, 'color', errbar_color{r}, 'linewidth', 1.1 );
            set( get(get(ebh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            % increase width of errorbars
            adjust_errorbars( ebh, plot_xlim );
        end
    end

    % add legend
    if (n==1)
        lh = legend( expt_strings );
        set(lh,'Location',legend_location);
    end

    %set(sph, 'XScale', 'log');
    box off;
    % set axes limits
    xlim( log10(plot_xlim) );
    ylim( plot_ylims{o} );
    % annotation and set axes
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

