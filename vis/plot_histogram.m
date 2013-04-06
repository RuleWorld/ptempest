function [] = plot_histogram( samples, cfg, figidx );
%PLOT_HISTOGRAM plot histograms of parameter samples
%
%  [] = plot_histogram( samples, cfg )
%  [] = plot_histogram( samples, cfg, figidx )
%  
%  where 'samples' is a (C x P x S) array of samples,
%  with C=number of chains, P=number of parameters, S=number of samples;
%  'cfg' is the configuration struct, and 'figidx' is an optional figure
%  index.
%
%  Histograms are shown for one parameter at a time. Hit any key to see
%  the next histogram.


% minimum param range to plot
width = 1.6;
param_range = [ cfg.param_location - width*cfg.param_scale; cfg.param_location + width*cfg.param_scale ];

nchains = size(samples,1);
nparams = size(samples,2);
nsamples = size(samples,3);

prompt = sprintf('goto [%d-%d], [n]ext,[p]rev,[q]uit? ', 1,nparams);

if exist('figidx') fh = figure(figidx);
else fh = figure; end
set( fh, 'Color', [1 1 1] );


res=30;
p=1;
while (p<=nparams)

    if (cfg.param_scale(p) == 0)
        p=p+1;
        continue;
    end

    maxpdf = 0;
    minlb = Inf;
    maxub = -Inf;
    for c = 1:nchains
        lb = min( permute( samples(c,p,:), [3 2 1]) );
        ub = max( permute( samples(c,p,:), [3 2 1]) );
        edges = linspace( lb, ub, res )';

        minlb = min([lb, minlb]);
        maxub = max([ub, maxub]);

        [N] = histc( permute( samples(c,p,:), [3 2 1]), edges(1:end-1) );
        pdf = N/sum((edges(2:end)-edges(1:end-1)).*N);
        maxpdf = max( [maxpdf; pdf ] );
    
        subplot(1,nchains,c);
        bar(edges(1:end-1),pdf,'histc');
        xlabel( sprintf('%s', cfg.param_names{p}), 'Interpreter', 'none' );
        ylabel( 'pdf' );
        title( sprintf('chain %d', c) );
    end

    for c = 1:nchains
        subplot(1,nchains,c);
        hold on;
        lb = min([param_range(1,p),minlb]);
        ub = max([param_range(2,p),maxub]);
        x = linspace( lb, ub, res*10 );
        fx = cfg.param_pdf{p}(x);
        plot( x, fx, '-r', 'linewidth', 2 );
        hold off;
        axis( [ lb, ub, 0, maxpdf] );
    end

    valid = 0;
    while (~valid)
        key = input( prompt, 's' );
        if strcmp(key,'') %next
            p = p+1;
            valid=1;
        elseif strcmp(key,'q') %quit
            p=nparams+1;
            valid=1;
        elseif strcmp(key,'p') %previous
            p = p-1;
            valid=1;
        elseif strcmp(key,'n') %next
            p = p+1;
            valid=1;
        else
            idx = str2num(key);
            if ( idx==round(idx) )
                if ( 1<=idx & idx<=nparams) %goto index
                    p=idx;
                    valid=1;
                else
                    fprintf(1,'index is out of range, try again!\n' ); 
                end
            else
                fprintf(1,'invalid input, try again!\n' ); 
            end
        end
    end    
    
end
