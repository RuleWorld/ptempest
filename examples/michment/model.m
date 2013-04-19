function [err, timepoints, species_out, observables_out ] = model( timepoints, species_init, parameters, suppress_plot )
%MODEL Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   'model' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the CVode library interfaced
%   to MATLAB via the MEX interface. Before running this script, the model
%   source in file model_cvode.c must be compiled (see that file for details).
%   MODEL returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = model( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   timepoints      : column vector of time points returned by integrator.
%   species_init    : row vector of 4 initial species populations.
%   parameters      : row vector of 5 model parameters.
%   suppress_plot   : 0 if a plot is desired (default), 1 if plot is suppressed.
%
%   Note: to specify default value for an input argument, pass the empty array.
%
%   OUTPUTS:
%   --------
%   err             : 0 if the integrator exits without error, non-zero otherwise.
%   timepoints      : a row vector of timepoints returned by the integrator.
%   species_out     : array of species population trajectories
%                        (columns correspond to species, rows correspond to time).
%   observables_out : array of observable trajectories
%                        (columns correspond to observables, rows correspond to time).
%
%   QUESTIONS about the BNG Mfile generator?  Email justinshogg@gmail.com



%% Process input arguments

% define any missing arguments
if ( nargin < 1 )
    timepoints = [];
end

if ( nargin < 2 )
    species_init = [];
end

if ( nargin < 3 )
    parameters = [];
end

if ( nargin < 4 )
    suppress_plot = 0;
end


% initialize outputs (to avoid error msgs if script terminates early
err = 0;
species_out     = [];
observables_out = [];


% setup default parameters, if necessary
if ( isempty(parameters) )
   parameters = [ 4.6052, 2.3026, -4.6052, 0.0, 0.0 ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  |  size(parameters,2) ~= 5  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 5].\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  |  size(species_init,2) ~= 4  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 4].\n' );
    err = 1;
    return;
end

% setup default timepoints, if necessary
if ( isempty(timepoints) )
   timepoints = linspace(0,1000,100+1)';
end
% check that timepoints has proper dimensions
if (  size(timepoints,1) < 2  |  size(timepoints,2) ~= 1  )
    fprintf( 1, 'Error: size of timepoints argument is invalid! Correct size = [t 1], t>1.\n' );
    err = 1;
    return;
end

% setup default suppress_plot, if necessary
if ( isempty(suppress_plot) )
   suppress_plot = 0;
end
% check that suppress_plot has proper dimensions
if ( size(suppress_plot,1) ~= 1  |  size(suppress_plot,2) ~= 1 )
    fprintf( 1, 'Error: suppress_plots argument should be a scalar!\n' );
    err = 1;
    return;
end

% define parameter labels (this is for the user's reference!)
param_labels = { 'S0', 'E0', 'log_k1', 'log_k2', 'log_k3' };



%% Integrate Network Model
try 
    % run simulation
    [err, species_out, observables_out] = model_cvode( timepoints, species_init, parameters );
catch
    fprintf( 1, 'Error: some problem integrating ODE network! (CVODE exitflag %d)\n', err );
    err = 1;
    return;
end



%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
    observable_labels = { 'S', 'E', 'ES', 'P' };

    % construct figure
    plot(timepoints,observables_out);
    title('model observables','fontSize',14,'Interpreter','none');
    axis([0 timepoints(end) 0 inf]);
    legend(observable_labels,'fontSize',10,'Interpreter','none');
    xlabel('time','fontSize',12,'Interpreter','none');
    ylabel('number','fontSize',12,'Interpreter','none');

end



%~~~~~~~~~~~~~~~~~~~~~%
% END of main script! %
%~~~~~~~~~~~~~~~~~~~~~%



% initialize species function
function [species_init] = initialize_species( params )

    species_init = zeros(1,4);
    species_init(1) = params(1);
    species_init(2) = params(2);
    species_init(3) = 0;
    species_init(4) = 0;

end


end
