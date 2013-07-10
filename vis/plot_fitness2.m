function [fitpdf,edges] = plot_fitness2( params1, params2, cfg1, cfg2, figidx )
%PLOT_FITNESS2 Get fitness values for the parameter ensemble
%
%  [fitness] = plot_fitness2( params1, params2, cfg1, cfg2 )
%  [fitness] = plot_fitness2( params1, params2, cfg1, cfg2, figidx )
%
%  where 'params1' and 'params2' are (S x P) arrays of parameters, 'cfg1' and
%  'cfg2' are the configuration structs, and 'figidx' is an optional figure index.
%
%  S = number of samples, P = number of parameters


% permute params if the samples are lie along the 3rd index
if (size(params1,1)==1 & size(params1,3)>1 )
    params1 = permute(params1,[3 2 1]);
    params2 = permute(params2,[3 2 1]);
end

% define xlabel
xlabel_strings = {'energy','-logprior','-loglike'};
% define ylabel
ylabel_string = 'pdf';



% dimensions
S = size(params1,1);
% array for fitness values
fitness1 = zeros(S,3);
fitness2 = zeros(S,3);
% remove time penalty
cfg.time_penalty = 0;

% loop over parameter samples
for s = 1:S

    fitness1(s,1) = cfg.energy_fcn( params1(s,:) );
    fitness1(s,2) = -cfg.logpdf_prior_fcn( params1(s,:) );
    fitness1(s,3) = fitness(s,1) - fitness(s,2);

    fitness2(s,1) = cfg.energy_fcn( params1(s,:) );
    fitness2(s,2) = -cfg.logpdf_prior_fcn( params1(s,:) );
    fitness2(s,3) = fitness(s,1) - fitness(s,2);

    % show progress
    fprintf(1,'.');
    if mod(s, 72)==0
        fprintf(1,'\n');
    end

end
fprintf(1,'\n');


if exist('figidx') fh = figure(figidx);
else fh = figure; end
set( fh, 'Color', [1 1 1] );

res=24;
for k = 1:3

    lb = min( fitness(:,k) );
    ub = max( fitness(:,k) );
    edges = linspace( lb, ub, res )';

    [N] = histc( fitness(:,k), edges(1:end-1) );
    pdf = N/sum((edges(2:end)-edges(1:end-1)).*N);

    subplot(1,3,k);
    bar(edges(1:end-1),pdf,'histc');
    xlabel( sprintf('%s', xlabel_strings{k}), 'Interpreter', 'none' );
    ylabel( ylabel_string );
    
end

% end plot_distr function
end

