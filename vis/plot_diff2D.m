function [] = plot_diff2D( samples1, samples2, parname1, parname2, cfg)
%PLOT_HIST2D plot difference between two 2D parameter distributions
%
%  [] = plot_diff2D( samples1, samples2, parname1, parname2, cfg )
%  
%  where 'samples1' is a (S1 x P) array of samples,
%  'samples2' is a (S2 x P) array of samples
%  with P=number of parameters, S=number of samples;
%  parname1 and parname2 are the names of the parameters to plot;
%  and 'cfg' is the configuration struct.

gridsize = 16;
levels = 16;

% permute params if the samples are lie along the 3rd index
if (size(samples1,1)==1 & size(samples1,3)>1 )
    samples1 = permute(samples1,[3 2 1]);
end
if (size(samples2,1)==1 & size(samples2,3)>1 )
    samples2 = permute(samples2,[3 2 1]);
end


% find parameter index
xpar = cfg.param_map.(parname1);
ypar = cfg.param_map.(parname2);

% find min and max
minx = min([samples1(:,xpar);samples2(:,xpar)]);
maxx = max([samples1(:,xpar);samples2(:,xpar)]);

miny = min([samples1(:,ypar);samples2(:,ypar)]);
maxy = max([samples1(:,ypar);samples2(:,ypar)]);

% generate bin edges
xedges = linspace(minx, maxx, gridsize);
yedges = linspace(miny, maxy, gridsize);

% get axis labels
xstring = sprintf( '%s (%s)', parname1, cfg.param_defs{xpar}.units );
ystring = sprintf( '%s (%s)', parname2, cfg.param_defs{ypar}.units );

% build 2D histograms
N1 = hist2( samples1(:,xpar), samples1(:,ypar), xedges, yedges);
N2 = hist2( samples2(:,xpar), samples2(:,ypar), xedges, yedges);

% normalize
dx = xedges(2) - xedges(1);
dy = yedges(2) - yedges(1);
pdf1 = N1/sum(sum( N1*dx*dy ));
pdf2 = N2/sum(sum( N2*dx*dy ));

% computer average maximum of pdf1 and pdf2
maxpdf = (max(max(pdf1)) + max(max(pdf2)))/2;

% compute normalized difference
Z = (pdf2-pdf1); %2*(pdf2 - pdf1)./(pdf2 + pdf1);
Z(isnan(Z)) = 0;

% find max absolute Z
maxabsz = max(max(abs(Z)));

% contourplot

subplot(2,2,1);
contourf(xedges, yedges, Z, linspace(-maxabsz,maxabsz,levels), 'linecolor','none');
% add color bar
colorbar;
% label graph
title( 'difference in 2D marginal distibrutions', 'fontsize',16 );
xlabel( xstring, 'fontsize',14, 'interpreter','none');
ylabel( ystring, 'fontsize',14, 'interpreter','none');

subplot(2,2,3);
contourf(xedges, yedges, pdf1, levels, 'linecolor','none');
% add color bar
colorbar;
% label graph
title( 'sample1, 2D marginal distibrutions', 'fontsize',16 );
xlabel( xstring, 'fontsize',14, 'interpreter','none');
ylabel( ystring, 'fontsize',14, 'interpreter','none');

subplot(2,2,4);
contourf(xedges, yedges, pdf2, levels, 'linecolor','none');
% add color bar
colorbar;
% label graph
title( 'sample2, 2D marginal distibrutions', 'fontsize',16 );
xlabel( xstring, 'fontsize',14, 'interpreter','none');
ylabel( ystring, 'fontsize',14, 'interpreter','none');

