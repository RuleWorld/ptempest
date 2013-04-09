function [err] = plot_pike( params, cfg, figidx )
%PLOT_PIKE Simulate and plot parameter scan over ligand and receptor concentration
%
%  Plot the distribution of saturation binding curves and the receptor dimer
%  fractions for the parameter set 'param':
%
%    [] = plot_pike( params, cfg )
%    [] = plot_pike( params, cfg, figidx )
%
%  where 'params' is a (S x P) or (1 x P x S) array of parameters,
%  'cfg' is the configuration struct, and 'figidx' is an optional figure index.
%
%  S = number of samples, P = number of parameters


% permute params if the samples are lie along the 3rd index
if (size(params,1)==1 & size(params,3)>1 )
    params = permute(params,[3 2 1]);
end


% set upper and lower quantiles for distribution fill
upperquant    = 0.978;
midupperquant = 0.842;
medquant      = 0.5;
midlowerquant = 0.158;
lowerquant    = 0.022;
% number of horizontal subplots
n_horiz_subplots = 3;
% legend location
legend_location = 'NorthWest';


% Find and sort EGF doses (M) and EGFR count (/cell) 
egf  = zeros(1,length(cfg.data));
egfr = zeros(1,length(cfg.data));
for d = 1 : length(cfg.data)
    egf(d) = cfg.data{d}.egf;
    egfr(d) = cfg.data{d}.egfr;
end
egf  = sort(unique(egf));
egfr = sort(unique(egfr));
nL = length(egf);
nR = length(egfr);


% get index of observables to display
obsv_freelig = [cfg.obsv_map.LigFree];
obsv_disp = [cfg.obsv_map.RecBound, cfg.obsv_map.RecDimer];


% dimensions
N = length(obsv_disp); % display observables
O = cfg.nobsv;         % observables
S = size(params,1);    % param samples


% simulate experiments
expt = {};
for r = 1:nR

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

            % convert free ligand fraction to concentration
            obsv(obsv_freelig) = obsv(obsv_freelig) * egf(l);

            % save observables
            expt{r}.obsv(1,:,s) = obsv;
    
        end

        % compute quantiles
        expt{r}.obsv_min(l,:)   = quantile( expt{r}.obsv(1,:,:), 0.0, 3);
        expt{r}.obsv_lower(l,:) = quantile( expt{r}.obsv(1,:,:), lowerquant, 3);
        expt{r}.obsv_med(l,:)   = quantile( expt{r}.obsv(1,:,:), medquant, 3);
        expt{r}.obsv_upper(l,:) = quantile( expt{r}.obsv(1,:,:), upperquant, 3);
        expt{r}.obsv_max(l,:)   = quantile( expt{r}.obsv(1,:,:), 1.0, 3);

        % get lab measurements
        expt{r}.mean(l,:)   = cfg.data{n}.mean(end,:);
        % convert free ligand fraction to concentration
        expt{r}.mean(l,obsv_freelig) = expt{r}.mean(l,obsv_freelig) * egf(l);

        expt{r}.stderr(l,:) = cfg.data{n}.stdev(end,:) ./ sqrt(cfg.data{n}.nsamples(end,:));
        % convert free ligand fraction to concentration
        expt{r}.stderr(l,obsv_freelig) = expt{r}.stderr(l,obsv_freelig) * egf(l);

        % show progress
        fprintf(1,'.');
        if mod((r-1)*nL + l, 72)==0
            fprintf(1,'\n');
        end

    end

end

% open figure
if exist('figidx') fh = figure(figidx);
else fh = figure; end
set( fh, 'Color', [1 1 1] );
set( fh,'PaperPositionMode','auto');
set( fh,'PaperOrientation','landscape');
hold off;


% style for errorbars
errbar_style = {'o','o','o'};
errbar_color = {[0 0 0.7], [0 0.7 0], [0.7 0 0]};
% style/color for lines
line_style = {'-b','-g','-r'};
% color for fill
fill_style = {'b','g','r'};


% xlabel
xlabel_string = sprintf( 'log10 [%s]', cfg.obsv_names{obsv_freelig} );


% define ylabels
ylabel_strings = {};
for o = obsv_disp
    ylabel_strings{o} = sprintf( '%s (%s)', cfg.obsv_names{o}, cfg.obsv_units{o} );
end


% experiment stings
expt_strings = {};
for r = 1:nR
    expt_strings{r} = sprintf('R/cell=%dk', round(egfr(r)/1000) );
end


% x axis limits
min_freelig = Inf;
max_freelig = -Inf;
for r = 1:nR
    if min_freelig > min( expt{r}.obsv_min(:,obsv_freelig) )
        min_freelig = min( expt{r}.obsv_min(:,obsv_freelig) );
    end
    if max_freelig < max( expt{r}.obsv_max(:,obsv_freelig) )
        max_freelig = max( expt{r}.obsv_max(:,obsv_freelig) );
    end
end
plot_xlim = [min_freelig*10/11, max_freelig*11/10];
% y axis limits
for o = 1:O
    if ( isfield(cfg,'plot_ylims') )
        plot_ylims{o} = cfg.plot_ylims{o};
    else
        plot_ylims{o} = [0 1.05];
    end
end


% plot curves
for n = 1:N
    
    % get observable index
    o = obsv_disp(n);

    sph = subplot( ceil(N/n_horiz_subplots), min(n_horiz_subplots,N), n);
    for r = 1:nR

        objh = fill( log10([expt{r}.obsv_med(:,obsv_freelig); flipud(expt{r}.obsv_med(:,obsv_freelig))]), ...
                     [expt{r}.obsv_lower(:,o); flipud(expt{r}.obsv_upper(:,o)) ], ...
                     fill_style{r}, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.7 );
        set( get(get(objh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        hold on;
        plot( log10(expt{r}.obsv_med(:,obsv_freelig)), expt{r}.obsv_med(:,o), line_style{r}, 'linewidth', 1.1 );
       
    end

    for r = 1:nR
        % plot errorbars (if any)
        if (n~=1) continue; end;
        if ( any(~isnan(expt{r}.mean(:,o))) )
            ebh = errorbar( log10(expt{r}.mean(:,obsv_freelig)), expt{r}.mean(:,o), ...
                            expt{r}.stderr(:,o), expt{r}.stderr(:,o), ...
                            errbar_style{r}, 'color', errbar_color{r}, 'linewidth', 1.1 );
            set( get(get(ebh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
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


