function [key_struct] = run_ptempest(BNGL_Model_Name,Experiments, N_Swaps, Job_Name)

%The purpose  of this matlab script is to point to ptempest, point to the
%model, and also point to the output location. The output location will
%contain all of the custom scripts needed in order  to succesffuly run
%ptemptest. 

%--------------------------------------------------------------------------
%                Preventing and Notifying about Common Errors
%--------------------------------------------------------------------------
%The BNGL model name must be a string and not a character vector. If character 
%vector is provided, it will be converted to a string. 
if ~isstring(BNGL_Model_Name) 
    if ischar(BNGL_Model_Name)
        BNGL_Model_Name = string(BNGL_Model_Name); 
    else 
        error("Model name must be a string or char");
    end 
end 

%Error if no experiments are defined: 
if ~exist('Experiments','var')
    error(sprintf("No experiments defined. Must create a struct object in the form of:\n<experiment_name>.<par_to_update> = <new value>"))
end 

%Number of swaps
if ~exist('N_Swaps', 'var')
    N_Swaps=200;
end

%Output Location 
Output_Loc = "Fits";
if ~exist('Job_Name', 'var')
    Job_Name= '';
end
     
%Finding the location where the run_ptempest function is stored. 
loc_of_pTempestWrapper = mfilename('fullpath');
[fPath, ~, ~] = fileparts(loc_of_pTempestWrapper); 

%pTempest is located in the directory above the runptempest file and the
%pTempest location should be absolute paths and NOT relative
pTempest_Loc = fullfile(fPath,"..");

%BNGL File Location 

BNGL_Function_Location = fullfile(fPath,"..","..","BioNetGen-2.5.1");

%General Model Information
BNGL_Model_Loc  = ".";
if ~exist('BNGL_Model_Name','var')
    BNGL_Model_Name = "Model"; 
end 

%CSV Files Needed to run pTemptest
Use_CSV_Files = false;
Parameter_CSV_Loc = "Parameters.csv";
Observables_CSV_Loc = "Observables.csv";
Data_CSV_Loc = "Data.csv";
Heuristic_CSV_Loc = "Heuristics.xlsx";

%Excel Document for Fit information: 
Fit_Information_Loc = "FitInformation.xlsx";

%Other Options: 
Sim_DT = 1;

%--------------------------------------------------------------------------
%Everything above should be just variable declarations and everythign below
%should be actions needed in order to run ptemptest. 
%--------------------------------------------------------------------------
key_struct = init_config_defaults();

key_struct.sim_dt = Sim_DT;
key_struct.nswaps = N_Swaps; 

%Key ptempest folders: 
key_struct.ptempest_loc = pTempest_Loc;
key_struct.core_loc     = fullfile(pTempest_Loc,"core");
key_struct.dist_loc     = fullfile(pTempest_Loc,"core","distr");
key_struct.vis_loc      = fullfile(pTempest_Loc,"vis");

%Model Locations
key_struct.bngl_model_loc = fullfile(pwd,BNGL_Model_Loc); 
key_struct.bngl_model_name = BNGL_Model_Name; 

%BNGL Location
key_struct.bngl_function_location = BNGL_Function_Location; 

%Orginal Key CSV file location
key_struct.use_csv_files = Use_CSV_Files;
if Use_CSV_Files
    key_struct.original_parameter_loc  = fullfile(pwd,Parameter_CSV_Loc);
    key_struct.original_observable_loc = fullfile(pwd,Observables_CSV_Loc);
    key_struct.original_data_loc       = fullfile(pwd,Data_CSV_Loc); 
    key_struct.original_heuristic_loc  = fullfile(pwd,Heuristic_CSV_Loc);
else 
    key_struct.orginal_fit_information_loc = fullfile(pwd,Fit_Information_Loc); 
end 

%Check to ensure all of the key  folders exists: 
check_ptempest_locations(key_struct); 

%Create Folder where results will be stored: 
key_struct.output_location = create_output_folders(Output_Loc,Job_Name); 

%Copy Files to final output job location 
key_struct = copy_files(key_struct);

%Get matlab files from BNGL model  
key_struct = compile_model(key_struct);

%Add Path to model so current and latest model is used: 
set_paths(key_struct);

%Get observable names, parameter names, and parameter values from the .m
%file 
key_struct = get_information_from_model(key_struct); 

%Create the param_defs array using orginal values specified in BNGL
key_struct = create_param_defs(key_struct);

%Update param_defs array to values/ranges as inputted in the Parameter_CSV_Loc
%file 
key_struct = update_param_defs(key_struct);
key_struct = init_parameter_defs(key_struct.param_defs,key_struct);


%Create the obsv_defs struct from model
key_struct = create_obs_defs(key_struct);
key_struct = init_observable_defs( key_struct.obsv_defs, key_struct );

%Loading experiments in Key Structures: 
key_struct.experiment = Experiments; 

%Loading Specific Names for Axis 
if exist("ExperimentNames")
    key_struct.experiment_plotting_names = ExperimentNames;
else 
    key_struct.experiment_plotting_names = struct(); 
end 

%Loading Specific Names for Axis 
if exist("AnalyteNamesOutputs")
    key_struct.analyte_plotting_names = AnalyteNamesOutputs;
else 
    key_struct.analyte_plotting_names = struct(); 
end 

%Loading Specific Time Information 
if exist("TimeInformation")
    key_struct.time_information_plotting = TimeInformation;
else 
    key_struct.time_information_plotting = struct(); 
end 

%Setting list of experiments to run. This is different from the set of
%experiments needed to fit because prediction simulations could be checked
%for convergence. 
key_struct.simulations_to_run = fieldnames(Experiments);

%Preparing the Data: 
key_struct = prepare_data(key_struct);

%Simulate Function
key_struct.simulate_fit_trajectories = @(params) simulate_fit_trajectories(params,key_struct);

key_struct.get_post_simulation_modifications = @(simulation_data) get_post_simulation_modifications(simulation_data,key_struct);

key_struct.get_simulated_data = @(Obsverable,Time,Experiments,SimulationData) get_simulated_data(Obsverable,Time,Experiments,SimulationData,key_struct);

key_struct.get_simulated_data_for_fitting = @(Obsverable,SimulationData) get_simulated_data_for_fitting(Obsverable,SimulationData,key_struct);

key_struct.set_paths = @() set_paths(key_struct);
%Energy Function: 

key_struct = setup_default_functions( key_struct );

key_struct.energy_fcn = @(params) energy_gaussian_custom_v2(params, key_struct);

%Running the Fit
startdir=pwd;
cd(key_struct.output_location)
save("key_struct.mat","key_struct");
parallel_tempering(key_struct);
cd(startdir);

end

%--------------------------------------------------------------------------
%Functions
%--------------------------------------------------------------------------
%Set the Paths: 
function [] = set_paths(StructObj)
    addpath(StructObj.bngl_output)
    addpath(StructObj.core_loc)
    addpath(StructObj.vis_loc)
    addpath(StructObj.dist_loc)
end 

%Intialize the key Structures: 
function [cfg] = init_config_defaults()
%INIT_CONFIG_DEFAULTS set default configuration options
%
%   [cfg] = init_config_defaults()

    cfg = struct();

    % parallel temperating options
    cfg.jobname = 'pt';                    % job name, for file input/output
    cfg.shuffle  = 1;                      % shuffle random number streams (seed by clock)
    cfg.parallel = 0;                      % parallel? true/false
    cfg.maxlabs  = 1;                      % maximum number of labs for parallel computing
    cfg.nchains  = 4;                      % number of chains
    cfg.nswaps = 2000;                      % number of chain swaps
    cfg.nsteps = 25;                       % number of steps between chain swaps
    cfg.display_status_interval = 10;       % How often to display info
    cfg.save_progress_interval = 100;     % How often to save progress 
    cfg.adapt_relstep_interval = 100;      % How often to adapt relative step size
    cfg.adapt_relstep_rate = 0.20;         % relative step size adaption rate
    cfg.optimal_step_acceptance = 0.24;    % optimal rate of step acceptance
    cfg.adapt_beta_interval = 250;         % how often to adapt temperature gradient
    cfg.adapt_beta_rate = 0.04;            % beta adaption rate
    cfg.optimal_swap_acceptance = 0.24;    % optimal rate of swap acceptance
    cfg.adapt_last = 2900;                 % last adaption step
    cfg.min_adaption_factor = 0.80;        % minimum adaption factor
    cfg.max_adaption_factor = 1.25;        % maximum adaption factor
    cfg.energy_init_max = 1000;            % maximum allowed energy for initialization
    cfg.max_init_steps = 500;              % maximum attempts at initialization
    cfg.max_beta = 1.0;                    % maximum chain beta (inverse of minimum chain temperature)
    cfg.beta_init = 0.666;                 % beta initialization parameter
    cfg.relstep_init = 0.01;               % relstep initialization parameter
    cfg.big_energy = 1e29;                 % big energy value: return if integrator fails
    cfg.timepenalty = 0;                   % penalty for long integration times
    
    % config trajectory visualization
    cfg.time_units = 's';                  % time units
    cfg.sim_tstart = 0;                    % simulation start time, s
    cfg.sim_tstop  = 20;                   % simulation stop time, s
    cfg.sim_dt     = 0.5;                  % time step for trajectory, s

end 

%Function: Check for key locations in ptempest 
function [] = check_ptempest_locations(StructObj)
    
    %Checks to ensure ptempest folder exists: 
    if ~exist(StructObj.ptempest_loc,'dir')
        error("Location of ptempest is  not found! Check pTempest_Loc variable!")
    else 
        alert =        ptempest_check(StructObj.core_loc,"core"); 
        alert =alert + ptempest_check(StructObj.dist_loc,fullfile("core","distr")); 
        alert =alert + ptempest_check(StructObj.vis_loc,"vis"); 
        if alert > 0 
            error("ptempest folder found but key folders are missing (check above)." + ... 
                   "Check ptempest download and consider re-downloading it from: https://github.com/RuleWorld/ptempest")
        end 
    end 
end 

%Function: To direct user to what key folders are missing in ptempest and
%          alert the user that the program could be corrupted 
function [alert] = ptempest_check(FolderLoc,FolderName)
    if ~exist(FolderLoc,'dir')
        alert = 1; 
        fprintf("Check the ptempest folder, the " + FolderName + " folder is missing\n")
    else 
        alert = 0;
    end 
end 

%Function: Create Output Folders 
function [output_loc] = create_output_folders(MainOutLoc,JobName)
    if ~exist(MainOutLoc,'dir')
        mkdir(MainOutLoc);
        fprintf("Created Output Directory: " + MainOutLoc +"\n")
    else 
        fprintf("'" + MainOutLoc+ "'" + " location already exists and will continue with the existing folder \n")
    end 
    
    %Getting the time and date formatted for the output folder name. In
    %order to reuse names of folders but not override existing folders
    %(jobs) that have the same name, a timestamp is added in the following
    %way: 
    % <Year>-<Month>-<Day>_<Hour>-<Minute>-<Second>_<JobName>
    time = datetime(); 
    
    date     = sprintf('%04d-%02d-%02d', time.Year, time.Month, time.Day);
    run_time = sprintf('%02d-%02d-%2d', time.Hour, time.Minute, round(time.Second));
    
    full_name = date +"_"+ run_time;
    if (JobName~="")
        full_name= full_name + '_' + JobName;
    end
    
    output_loc = fullfile(MainOutLoc,full_name);
    if ~exist(output_loc,'dir')
        mkdir(output_loc);
    else 
        error("Final output folder exists. Check create_output_folders function")
    end 
    
    output_loc = fullfile(pwd, output_loc);
end 


%Function: Copy files needed to run the fit in its current condition 
function [StructObj] = copy_files(StructObj) 
    loc_of_bngl_model = fullfile(StructObj.bngl_model_loc, ... 
                                 StructObj.bngl_model_name) + ".bngl"; 
                    
    if StructObj.use_csv_files
        data_csv       = StructObj.original_data_loc;
        observable_csv = StructObj.original_observable_loc;
        parameter_csv  = StructObj.original_parameter_loc;
        heuristics_csv = StructObj.original_heuristic_loc;
        files_2_copy = [data_csv observable_csv parameter_csv ... 
                        heuristics_csv loc_of_bngl_model];

        names = ["data_loc","observable_loc","parameter_loc",...
                 "heuristic_loc","model_loc"];
    else 
        orginal_fit_information = StructObj.orginal_fit_information_loc;
        files_2_copy = [orginal_fit_information loc_of_bngl_model];
        names = ["fit_information_loc","model_loc"];
    end 
        
                
    for i_files = 1:length(files_2_copy)
        copy_file_i = files_2_copy{i_files}; 
        
        %From the file path, obtain the name of the file to be copied 
        [~,name,ext] = fileparts(copy_file_i);
        final_file_name = strcat(name,ext);
        StructObj.(names(i_files)) = fullfile(StructObj.output_location,final_file_name);
        [status,~] = copyfile(copy_file_i,StructObj.(names(i_files)));
        if ~status
            error("The following file was not found: " + copy_file_i)
        end 
    end 
end 

%Function: Run BNGL model and compile if necessary to get .m file 
function [StructObj] = compile_model(StructObj) 
    full_function = fullfile(StructObj.bngl_function_location,"BNG2.pl");
    model_run = fullfile(StructObj.output_location,strcat(StructObj.bngl_model_name,".bngl"));
    
    %Run BNGL script to generate matlab files
    [~,status] = perl(full_function,model_run); 
    
    if status 
        error("BNGL model failed to run and make .m file, check model file"); 
    end 
    
    %Make Output Directo and Move Files to that location 
    out_location = fullfile(StructObj.output_location,"BNGL_Output");
    mkdir(out_location);
    
    %Assume the code needs compiled unless _cvode.c is not found 
    StructObj.compile_code = 1; 
    files_2_move = [".cdat" ".gdat" ".m" ".net" "_cvode.c"];
    
    for i_files_2_move = 1:length(files_2_move) 
        file_n = StructObj.bngl_model_name + files_2_move{i_files_2_move};
        status = movefile(file_n, out_location);
        if ~status && files_2_move{i_files_2_move} == "_cvode.c"
            StructObj.compile_code = 0; 
        end 
    end 
    
    StructObj.bngl_output = out_location;
    
    %Compile Model if writeMexfile was used: 
    if StructObj.compile_code
        compile_mex_model(StructObj);
    end 
end 

%Function: Compile the mex model 
function [] = compile_mex_model(StructObj)
    
    loc_of_includeFile = fullfile(StructObj.bngl_function_location, "include");
    
    loc_of_lib_file = fullfile(StructObj.bngl_function_location,"lib");
    
    cvode_loc = fullfile(StructObj.bngl_output,strcat(StructObj.bngl_model_name,"_cvode.c"));

    CompileCode = ...
         sprintf('mex %s -I%s -L%s -lsundials_cvode -lsundials_nvecserial',cvode_loc,loc_of_includeFile,loc_of_lib_file);
    eval(CompileCode);

    movefile(StructObj.bngl_model_name + "_cvode.mex*",StructObj.bngl_output);
end 

%Function: Get model information: 
function [StructObj] = get_information_from_model(StructObj)
    model = fullfile(StructObj.bngl_output,StructObj.bngl_model_name);
    % Get parameter names from param_labels string definition
    [~, result]= system(sprintf('grep "param_labels = { " %s.m', model));
    eval(result);
    StructObj.param_labels = param_labels;
    
    % Get default parameter values from model 
    [~, result]= system(sprintf('grep "parameters = \\[ " %s.m', model));
    eval(result);
    StructObj.parameters = parameters;
    
    [~, result]= system(sprintf('grep "observable_labels = { " %s.m', model));
    eval(result);
    StructObj.observable_labels = observable_labels;
end 

%Function: Create the param_defs cell array assuming all inputs are point 
%          distributions set at the same value as in BNGL
function [StructObj] = create_param_defs(StructObj)
    n_parameters = length(StructObj.parameters);
    param_defs = cell(1,n_parameters);
    
    param_index = struct();
    
    for ith_parameter = 1:n_parameters
        temp_name = StructObj.param_labels{ith_parameter};
        temp_param_def = struct();
        temp_param_def.name  = temp_name;
        temp_param_def.prior = "point"; 
        temp_param_def.value = StructObj.parameters(ith_parameter); 
        param_defs{ith_parameter} = temp_param_def;
        param_index.(temp_name) = ith_parameter;
    end 
    
    StructObj.original_param_defs = param_defs;
    StructObj.param_index = param_index;
end


%Function: Create the param_defs file by updating the pars as specified by
%          an external excel sheet 
function [StructObj] = update_param_defs(StructObj)
    param_defs = StructObj.original_param_defs; 
    
    %Load parameter changes file and ensure there are no duplicate entries 
    if StructObj.use_csv_files
        parameter_changes = readtable(StructObj.parameter_loc,"PreserveVariableNames",true);
    else
        parameter_changes =readtable(StructObj.fit_information_loc,"Sheet","Parameters","PreserveVariableNames",true);
    end 
    
    [n_updated_pars,~] = size(parameter_changes);
    
    [n_unique_pars,~] = size(unique(parameter_changes.Name));
    
    if n_updated_pars ~= n_unique_pars
        error("Check the parameters file, duplicate name found")
    end 
    
    %Update Parameters as specified in the parameter_changes file 
    for ith_update = 1:n_updated_pars
        par_to_update = parameter_changes(ith_update,:); 
        
        name_of_par_to_update = par_to_update.Name{1};
        
        try
            index_of_update = StructObj.param_index.(name_of_par_to_update); 
        catch 
            error("The following is not a parameter in the model: "+name_of_par_to_update)
        end 
        
        previous_param_def = param_defs{index_of_update};
        
        if ~strcmp(previous_param_def.name,name_of_par_to_update)
            error("Something is wrong with finding the correct index of the parameter")
        end 
        
        temp_param_def = struct(); 
        
        temp_param_def.name = par_to_update.Name{1};
        temp_param_def.units = par_to_update.units{1};
        temp_param_def.prior = par_to_update.prior{1}; 
        temp_param_def.(par_to_update.arg1{1}) = par_to_update.value1; 
        
        temp_param_def = test_if_na(temp_param_def,par_to_update.arg2,par_to_update.value2);
        temp_param_def = test_if_na(temp_param_def,par_to_update.arg3,par_to_update.value3);
        temp_param_def = test_if_na(temp_param_def,par_to_update.arg4,par_to_update.value4);
        
        param_defs{index_of_update} = temp_param_def;
    end 
    
    %Sub-Function: Not all priors have 4 arguments, therefore this function
    %is set to ignore those arguments if they do not exists
    function [temp_param_def] = test_if_na(temp_param_def,arg,value)
        arg_1_remove_spaces = strrep(arg{1}," ","");
        if  arg_1_remove_spaces ~='N/A'
            temp_param_def.(arg_1_remove_spaces) = value;
        end 
    end 

    StructObj.parameter_changes = parameter_changes;
    StructObj.param_defs = param_defs;
end

%Function: Create the observable labels and update units if the observable
%exists in the observable csv. Honestly, this step could be removed
%because I dont see the point of it. 
function [StructObj] = create_obs_defs(StructObj) 
    n_obs = length(StructObj.observable_labels);
    
    obsv_defs = cell(1,n_obs); 
    obsv_index = struct();
    
    for ith_obs = 1:n_obs
        temp_obs_struct = struct(); 
        temp_obs_struct.name = StructObj.observable_labels{ith_obs};
        temp_obs_struct.units = "N/A"; 
        temp_obs_struct.display = 1; 
        temp_obs_struct.minplot = 0;
        temp_obs_struct.maxplot = 1;
        
        obsv_defs{ith_obs} = temp_obs_struct;
        obsv_index.(StructObj.observable_labels{ith_obs}) = ith_obs;
    end 
    
    StructObj.obsv_defs = obsv_defs;  
    StructObj.n_obs = n_obs;
    StructObj.obsv_index = obsv_index;
end 

%Function: Prepare the data. 
function [StructObj] = prepare_data(StructObj)
    %Loading Data
    if StructObj.use_csv_files
        raw_data = readtable(StructObj.data_loc,"PreserveVariableNames",true);
    else
        raw_data = readtable(StructObj.fit_information_loc,"Sheet","Data","PreserveVariableNames",true);
    end 
    
    %Data Extraction
    analyte_names     = unique(raw_data.analyte);
    times_unique      = unique(raw_data.time);
    
    %Data by the numbers 
    n_fitted_analytes = length(analyte_names);
    
    %Determining the Simulation Time for the experiments: 
    StructObj.simulation_time = 0:StructObj.sim_dt:max(times_unique);
    
    %Extract Data Conversions if Necessary:
    post_simulation_analytes= struct();
    data = struct();
    
    for ith_analyte = 1:n_fitted_analytes
        
        %Get analyte that needs extracted 
        analyte_i = analyte_names{ith_analyte};
        
        %If the analyte has a triple underscore, that means  the analyte
        %requires a post simulation modification. The type and observable
        %are recorded. The type must correspond to a function that will
        %calculate the post simulation conversion parameters. 
        [split_names,matches] = split(analyte_i,"___");
        if ~isempty(matches)
            post_simulation_analytes.(analyte_i).type = split_names{2};
            post_simulation_analytes.(analyte_i).obsv = split_names{1};
        end 
        
        %Extracting the data that corresponds to analyte_i from the raw
        %config file. 
        reduced_data_analyte_i          = raw_data(strcmp(analyte_i,raw_data.analyte),:);
        
        %Creating a matrix where the x-axis represets time and each column
        %is a different experiment. 
        %Important Note: Assumes each time/experiment only contains a
        %                single entry. 
        value      = unstack(reduced_data_analyte_i,'value','experiment',"GroupingVariables",'time');
        sigma      = unstack(reduced_data_analyte_i,'sigma','experiment',"GroupingVariables",'time');
        weight     = unstack(reduced_data_analyte_i,'weight','experiment',"GroupingVariables",'time');
        
        
        %The first column is the time column and it will be recorded in a
        %separate struct object. The value, sigma, and weight will be
        %stored as an array object for easy computation of the energy
        %function. 
        data.(analyte_i).value      = table2array(value(:,2:end)); 
        data.(analyte_i).sigma      = table2array(sigma(:,2:end)); 
        data.(analyte_i).weight     = table2array(weight(:,2:end)); 
        data.(analyte_i).experiment = value.Properties.VariableNames(2:end);
        data.(analyte_i).time       = value.time;
        
        %Storing the Indices of the relavant timepoints 
        data.(analyte_i).time_index       = get_index(StructObj.simulation_time,value.time);
        data.(analyte_i).experiment_index = get_index(StructObj.simulations_to_run,data.(analyte_i).experiment);
    end 
    
    %Parameters to Save
    data.raw_data = raw_data; 
    data.analyte_names = analyte_names; 
    data.n_fitted_analytes = n_fitted_analytes; 
    
    
    StructObj.data = data;
    StructObj.post_simulation_analytes = post_simulation_analytes;
    
    %Load Heuristics: 
    if StructObj.use_csv_files
        StructObj.heuristic_table = readtable(StructObj.heuristic_loc,'PreserveVariableNames',true);
    else
        StructObj.heuristic_table = readtable(StructObj.fit_information_loc,"Sheet","Heuristics","PreserveVariableNames",true);
    end 
end

function [obsverable_out,err] = simulate_fit_trajectories(Parameters,StructObj)
    obsverable_out = zeros(length(StructObj.simulation_time),StructObj.n_obs,length(StructObj.simulations_to_run));
    
    for ith_experiment = 1:length(StructObj.simulations_to_run)
        experiment_i = StructObj.simulations_to_run{ith_experiment};
        [err,~,obsv] = simulate_model(StructObj.simulation_time',[],Parameters,StructObj.param_index,StructObj.experiment.(experiment_i),StructObj.bngl_model_name);
        
        if (err)
            obsverable_out = [];
            return;
        end
        
        obsverable_out(:,:,ith_experiment) = obsv;
    end 
    
end 


function [err, sp, obsv] = simulate_model( t, init, params,param_index,exp,model_name)

    pars_to_update = fieldnames(exp);

    for ith_par_update = 1:length(pars_to_update)
        par_name = pars_to_update{ith_par_update};
        params(param_index.(par_name)) = exp.(par_name);
    end 

    % run simulation
    model_2_run = str2func(model_name);
    [err, ~, sp, obsv] = model_2_run( t, init, params, 1 );
    if (err)
        sp = [];
        obsv = [];
        return;
    end

end 



function [PSMD] = get_post_simulation_modifications(SimulationData,StructObj)
    %Information Needed from the StructObj 
    post_simulation_analytes = StructObj.post_simulation_analytes;
    simulations = StructObj.simulations_to_run;
    simulation_time = StructObj.simulation_time;
    
    %Extracted Information 
    names_of_psa             = fieldnames(post_simulation_analytes);
    n_psa                    = length(names_of_psa) ;
    n_simulations            = length(simulations);
    n_time                   = length(simulation_time);
    
    
    PSMD = struct();
    %Looping through each post simulation analyte to calculate the
    %conversion factors 
    for ith_psa = 1:n_psa
        psa_i = names_of_psa{ith_psa};
        
        type = post_simulation_analytes.(psa_i).type;
        obsv = post_simulation_analytes.(psa_i).obsv;
        
        times_to_extract = StructObj.data.(psa_i).time_index;
        experiments_to_extract= StructObj.data.(psa_i).experiment_index; 
        obsv_index = StructObj.obsv_index.(obsv);
        exp_values = StructObj.data.(psa_i).value; 
        
        obsv_for_conversion = squeeze(SimulationData(times_to_extract,obsv_index,experiments_to_extract));
        
        conversion_function = str2func(type);
        
        parameters = conversion_function(obsv_for_conversion,exp_values);
        
        
        x = squeeze(SimulationData(:,StructObj.obsv_index.(obsv),:));
        x_post_modification = polyval(parameters,x); 
        
        PSMD.(psa_i).values = x_post_modification; 
        PSMD.(psa_i).parameters =parameters; 
        
        
    end 
    

end 
%Function: Get the data for fitting the model with the results. The
%function below will have the same output as get_simulated_data except it
%will not need to find the index locations of the times and experiments
%because it will already be known. This is in hopes to speed up the
%parallel tempering process and avoid unnecessary calculations. 
function [SimulationOut] = get_simulated_data_for_fitting(Obsverable,SimulationData,StructObj)
    %Finding the index that are associated in the data. The simulated
    %data will be extracted in the same order as the inputs of Time and experiments. 
    if length(size(SimulationData)) == 2
        index_of_obsv = 'N/A';
    else 
        index_of_obsv = StructObj.obsv_index.(Obsverable);
    end 
    index_of_experiments = StructObj.data.(Obsverable).experiment_index;
    index_of_time = StructObj.data.(Obsverable).time_index;
    
    if isstr(index_of_obsv)
        reduced_sim_data = SimulationData(index_of_time,index_of_experiments);
    else
        reduced_sim_data = SimulationData(index_of_time,index_of_obsv,index_of_experiments);
    end 

    SimulationOut = squeeze(reduced_sim_data);
end 


%Function: Accepts an observable label, a list of times, a list of
%experiments, and simulation data and returns a matrix where each row
%represents a time and column represents an experiment. The order is in the
%same order as the inputs 
function [SimulationOut] = get_simulated_data(Obsverable,Time,Experiments,SimulationData,StructObj)
    %Finding the index that are associated in the data. The simulated
    %data will be extracted in the same order as the inputs of Time and experiments. 
    if length(size(SimulationData)) == 2
        index_of_obsv = 'N/A';
    else 
        index_of_obsv = get_index(StructObj.observable_labels,Obsverable);
    end 
    index_of_experiments = get_index(StructObj.simulations_to_run,Experiments);
    index_of_time = get_index(StructObj.simulation_time,Time);
    
    if isstr(index_of_obsv)
        reduced_sim_data = SimulationData(index_of_time,index_of_experiments);
    else
        reduced_sim_data = SimulationData(index_of_time,index_of_obsv,index_of_experiments);
    end 

    SimulationOut = squeeze(reduced_sim_data);
end 

%Get the index locations of the orginal list of the specificed new list
%items. Index locations should be returned in the same order as the NewList
function [index_out] = get_index(OrginalList,NewList)
    [x,ia,ib] = intersect(OrginalList,NewList);
    
    %If the length of x does not match the new list, and it is not a
    %character, then the intersect failed ot find all of the values.
    %Sometimes a single 'char' list will be used as an input instead of a
    %cell array which causes the length of the NewList value to be function
    %of the number of characters which is wrong. 
    if length(x) ~= length(NewList) && ~ischar(NewList)
        index_out = [];
        for ith_time = 1:length(NewList)
            [~,closest_index] = min(abs(NewList(ith_time)-OrginalList));
            index_out = [index_out closest_index];
            
        end 
    else 
        %ia and ib are sorted but I want the index to return in the same
        %order as inputted by the New List, therefore I did the following to
        %make that happen. 
        t  = sortrows([ib ia ]);
        index_out = t(:,2);
    end 
end    

%Post Modification conversion when the conversion factor is relative to max
function [conversion] = relative_to_max(SimData,~)
        conversion= [1/max(SimData,[],'all') 0];
end 

%Post Modification conversion when the conversion factor is the best linear
%fit. 
function [conversion] = mle(SimData,ExpData)
    [rows,columns] = size(SimData);
    
    x = reshape(SimData,[1,rows*columns]);
    y = reshape(ExpData,[1,rows*columns]);
    
    conversion = polyfit(x,y,1);
end 

%Function: Energy Function Calculation: 
function [energy] = energy_gaussian_custom_v2( params, cfg )

    %ENERGY_GAUSSIAN Calculate parameter energy for gaussian likelihood 
    %model
    %  [energy] = energy_generic( params, cfg )

    %======================================================================
    % check prior probability on parameters
    %======================================================================
    logprior = cfg.logpdf_prior_fcn(params);
    if isinf(logprior)
        energy = cfg.big_energy;
        return;
    end

    energy = -logprior;
    
    %======================================================================
    % Starting the Simulation Timer
    %======================================================================
    % start integration timer
    simtimer = tic;
    
    %======================================================================
    % Simulate the Experiments
    %======================================================================
    
    [simulation_data,err] = cfg.simulate_fit_trajectories(params);
    if (err)
        energy = cfg.big_energy;
        return;
    end

    
    %======================================================================
    % Ending the Simulation Timer
    %======================================================================
    % penalize for slow integrations
    dt = toc(simtimer);
    
    %======================================================================
    % Post Process the Data
    %======================================================================
    post_sim_modification = cfg.get_post_simulation_modifications(simulation_data);
    
    %======================================================================
    % Calculate Energy Function
    %======================================================================
    for ith_analytes = 1:cfg.data.n_fitted_analytes
        analyte_i = cfg.data.analyte_names{ith_analytes};
        
        data_to_fit = cfg.data.(analyte_i);
        
        if isfield(cfg.post_simulation_analytes,analyte_i)
            sim_data = cfg.get_simulated_data_for_fitting(analyte_i,post_sim_modification.(analyte_i).values);
        else 
            sim_data = cfg.get_simulated_data_for_fitting(analyte_i,simulation_data);
        end 
        
        loglike = nansum(nansum( -data_to_fit.weight .* (sim_data - data_to_fit.value).^2 ./ (2*(data_to_fit.sigma).^2) ));
        energy = energy - loglike;
    end 
    
    %======================================================================
    % Calculate Heuristic Energy Function
    %======================================================================
    energy = energy + GetHeuristicTotals(cfg,simulation_data);
end 




%==========================================================================
%Functions Needed to Calculate Heuristics 
%==========================================================================
function [TotalHeuristicEnergy] = GetHeuristicTotals(StructObj,Simulation)
    TotalHeuristicEnergy = 0; 
    heuristic_table = StructObj.heuristic_table;

    [n_heuristics,~] = size(heuristic_table);

    for ith_heuristic = 1:n_heuristics
        heuristic_i = heuristic_table(ith_heuristic,:); 

        value_i = getValue(Simulation,heuristic_i,StructObj);

        energy = calcEnergy(value_i,heuristic_i);
        TotalHeuristicEnergy = TotalHeuristicEnergy + energy;
    end 

end 





function [Energy] = calcEnergy(Value,Heuristic)
    max_value = Heuristic.("Max Value");
    min_value = Heuristic.("Min Value"); 
    Weight = Heuristic.Weight;
    if ( min_value< Value  ) & ( Value < max_value ) 
        Energy = 0; 
    elseif min_value> Value
        Energy = (min_value - Value)*Weight;
    else 
        Energy = (Value - max_value)*Weight;
    end 
end 

function [ValueOut] = getValue(Simulation,heuristic_i,key_struct)
    experiment = heuristic_i.Experiment;
    if ~isnan(heuristic_i.Time)
        time = heuristic_i.Time;
    else 
        time = key_struct.simulation_time;
    end
    
    type = str2func(heuristic_i.Type{1});
    
    [numerator,denominator] = assess_observable(heuristic_i.Observable);
    
    

    numerator_values = key_struct.get_simulated_data(numerator,time,experiment,Simulation);
    denom_values = key_struct.get_simulated_data(denominator,time,experiment,Simulation);
    
    
    constants_numer = getConstants(setdiff(numerator,key_struct.observable_labels));
    constants_denom = getConstants(setdiff(denominator,key_struct.observable_labels));
    
    
    CalcValue = (constants_numer+numerator_values)./(constants_denom+sum(denom_values,2));
    
    ValueOut = type(CalcValue);


end 

function [out] = getConstants(InputCell)
    out = 0; 
    if ~isempty(InputCell)
        for ith_constant = 1:length(InputCell)
            out = out + str2num(InputCell{ith_constant});
        end 
    end 
end 

function [numerator,denominator] = assess_observable(Observable)
    if iscell(Observable)
        Observable = Observable{1};
    end 
    
    %Remove Empty Spaces 
    Observable = replace(Observable," ","");
    
    %Check for Division: 
    [numerator,denominator] = check_for_division(Observable);
    
    %Remove brackets: 
    numerator = remove_brackets(numerator);
    denominator = remove_brackets(denominator);
    
    %Split Summations 
    numerator = split_sums(numerator);
    denominator = split_sums(denominator);
    
    
end 

function [numerator,denominator] = check_for_division(input_observable)
    input_observable = split(input_observable,"/");
    if length(input_observable)==2 
        numerator = input_observable{1};
        denominator = input_observable{2};
    else 
        numerator = input_observable{1};
        denominator = '1';
    end
end 

function [output] = remove_brackets(input_equation) 

    output = replace(input_equation,"(","");
    output = replace(output,")","");
end 

function [output] = split_sums(input_equation) 
    output = split(input_equation,"+");
end 




