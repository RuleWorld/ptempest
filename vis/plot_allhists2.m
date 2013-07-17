function [] = plot_allhists2( samples1, samples2, cfg, fighandle, figname );
%PLOT_ALLHISTS2 plot all histograms of parameter samples
%
%  [] = plot_allhists2( samples1, samples2, cfg )
%  [] = plot_allhists2( samples1, samples2, cfg, fighandle )
%  [] = plot_allhists2( samples1, samples2, cfg, fighandle, figname )
%  
%  where 'samples1' and 'samples2' are (S x P) arrays of samples,
%  with C=number of chains, P=number of parameters, S=number of samples;
%  'cfg' is the configuration struct, and 'fighandle' is an optional figure
%  handle.


if and( size(samples1,1) == 1, size(samples1,3) > 1)
    samples1 = permute(samples1,[3 2 1]);
    samples2 = permute(samples2,[3 2 1]);
end

% minimum param range to plot
width = 1.6;
param_range = [ cfg.param_location - width*cfg.param_scale; cfg.param_location + width*cfg.param_scale ];

nsamples = size(samples1,1);
nparams = size(samples1,2);

if exist('fighandle')
    fh = fighandle;
else
    fh = gcf;
end
set( fh, 'Color', [1 1 1] );


plotparams = find(cfg.param_scale ~= 0);
Nplotparams = length(plotparams);

ncols = 5;
nrows = ceil(Nplotparams/ncols);

res=16;
for n = 1:Nplotparams

    p = plotparams(n);

    maxpdf = 0;

    lb = min( [samples1(:,p); samples2(:,p)] );
    ub = max( [samples1(:,p); samples2(:,p)] );
    edges = linspace( lb, ub, res )';

    % get histogram
    [N1] = histc( samples1(:,p), edges );
    [N2] = histc( samples2(:,p), edges );
    % get singleton in the last bin and put in previous bin
    N1(end-1) = N1(end-1)+N1(end);
    N2(end-1) = N2(end-1)+N2(end);
    % now shorten N by 1.
    N1 = N1(1:end-1);
    N2 = N2(1:end-1);

    % normalize histogram to get approx pdf
    pdf1 = N1/sum((edges(2:end)-edges(1:end-1)).*N1);
    pdf2 = N2/sum((edges(2:end)-edges(1:end-1)).*N2);
    maxpdf = max([pdf1 pdf2]);

    % setup subplot
    sph = subplot(nrows,ncols,n);
    %set( sph, 'Position', [0.20 0.20 0.75 0.75] );
    set( sph, 'Box', 'off');
    set( sph, 'fontsize', 8);

    lb = min([param_range(1,p),lb]);
    ub = max([param_range(2,p),ub]);
    x = linspace( lb, ub, res*10 )';
    fx = cfg.param_pdf{p}(x);
    maxpdf = max([fx; pdf1; pdf2]);
    objh = plot( x, fx, 'linestyle', '-', 'color', 0.8*[237/255,187/255,78/255], 'linewidth', 1.5 );
    set( get(get(objh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on;

    bh1 = bar(edges(1:end-1),pdf1,'histc');
    set(bh1,'facealpha',0.6);
    set(bh1,'facecolor',[0.8 0 0]);
    set(bh1,'edgecolor',[0.5 0 0]);

    bh2 = bar(edges(1:end-1),pdf2,'histc');
    set(bh2,'facealpha',0.6);
    set(bh2,'facecolor',[0 0 0.8]);
    set(bh2,'edgecolor',[0 0 0.5]);

    xlh = xlabel( sprintf('%s', cfg.param_names{p}), 'Interpreter', 'none','fontsize',10 );
    set(xlh, 'Position', [(ub+lb)/2,-0.18*maxpdf]);

    xticks = [lb,ub];
    mult = 1;
    while ceil(xticks(1)) >= floor(xticks(2))
        xticks = xticks*10;
        mult = mult*10;
    end
    xticks = [ceil(xticks(1))/mult,floor(xticks(2))/mult];

    if any( abs(xticks - round(xticks)) > 1e-8 )
        xticklabels = sprintf('%.1f|',xticks);
    else
        xticklabels = sprintf('%.0f|',xticks);
    end

    yticks = [0,maxpdf];
    mult = 1;
    while ceil(yticks(1)) >= floor(yticks(2))
        yticks = yticks*10;
        mult = mult*10;
    end
    yticks = [0,floor(yticks(2))/mult];

    if any( abs(yticks - round(yticks)) > 1e-8 )
        yticklabels = sprintf('%.1f|',yticks);
    else
        yticklabels = sprintf('%.0f|',yticks);
    end

    set(sph, ...
        'XTick',      xticks, ...
        'XTickLabel', xticklabels, ...
        'YTick',      yticks, ...
        'YTickLabel', yticklabels, ...
        'Box',        'off', ...
        'TickDir',    'out', ...
        'TickLength', [.04 .04], ...
        'XMinorTick', 'off', ...
        'YMinorTick', 'off', ...
        'XGrid',      'off', ...
        'YGrid',      'off' ...
    );

    axis([ lb, ub, 0, maxpdf]);
    hold off;

end

outfile = '';
if exist('figname')
    outfile = sprintf('%s_allhists', figname);
else
    outfile = 'allhists';
end

% finalize size of figure
set( fh, 'WindowStyle', 'normal' );
set( fh, 'Units', 'inches' );
set( fh, 'Position', [1 1 10 6] );

set( fh, 'PaperPositionMode','auto');
set( fh, 'PaperSize',[10 6]);
%print(fh, outfile, '-dpdf')


