function [] = plot_hist2D( samples, parname1, parname2, cfg)
%PLOT_HIST2D plot estimate of 2D parameter marginal distribution
%
%  [] = plot_hist2D( samples, parname1, parname2, cfg )
%  
%  where 'samples' is a (S x P) array of samples,
%  with P=number of parameters, S=number of samples;
%  parname1 and parname2 are the names of the parameters to plot;
%  and 'cfg' is the configuration struct.

gridsize = 15;
levels = 15;

% permute params if the samples are lie along the 3rd index
if (size(samples,1)==1 & size(samples,3)>1 )
    samples = permute(samples,[3 2 1]);
end

% find parameter index
xpar = cfg.param_map.(parname1);
ypar = cfg.param_map.(parname2);

% find min and max
minx = min(samples(:,xpar));
maxx = max(samples(:,xpar));

miny = min(samples(:,ypar));
maxy = max(samples(:,ypar));

% generate bin edges
xedges = linspace(minx, maxx, gridsize);
yedges = linspace(miny, maxy, gridsize);

% get axis labels
xstring = sprintf( '%s (%s)', parname1, cfg.param_defs{xpar}.units );
ystring = sprintf( '%s (%s)', parname2, cfg.param_defs{ypar}.units );

% build 2D histogram
N = hist2( samples(:,xpar), samples(:,ypar), xedges, yedges);

% normalize
dx = xedges(2) - xedges(1);
dy = yedges(2) - yedges(1);
pdf = N/sum(sum( N*dx*dy ));

% contourplot
contourf(xedges, yedges, pdf, levels);
% add color bar
colorbar;
% label graph
title( '2D marginal distibrution, sample estimate', 'fontsize',16 );
xlabel( xstring, 'fontsize',14, 'interpreter','none');
ylabel( ystring, 'fontsize',14, 'interpreter','none');

