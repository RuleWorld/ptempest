classdef plotting_tool < handle
   properties
      loc_key_struct 
      num_samples    
      converged_start
      
      %Temporary Visable Parameters:
      key_struct
      last_swap
      samples_chosen
      fitted_simulated_data
   end
   
   properties(Hidden)
      params_chain
      energy_chain
      current_converged_start
      %Parameters useful to be visable while developing:
      
   end 
   
   methods
       function [obj] = sample_patients(obj,ConvergedStart)
           if ~exist('ConvergedStart','var')
               ConvergedStart = obj.converged_start; 
           else 
               obj.converged_start = ConvergedStart;
           end
           
           obj.samples_chosen = datasample(ConvergedStart:obj.last_swap,...
                                 obj.num_samples,"Replace",false);
           obj.current_converged_start = ConvergedStart;
       end 
       function [] = plot_all_energy_chains(obj,Start,RollingAvg,Stop)
            if ~exist('Start','var')
                Start = 1; 
            end

            if ~exist('Stop','var')
                Stop = obj.last_swap; 
            end

            if ~exist('RollingAvg','var')
                RollingAvg = 1; 
            end 
            energy_chain_to_plot = obj.energy_chain;
            plot_energy_chain(obj,energy_chain_to_plot,Start,RollingAvg,Stop)
       end 
       function [] = plot_top_energy_chain(obj,Start,RollingAvg,Stop)
            if ~exist('Start','var')
                Start = 1; 
            end

            if ~exist('Stop','var')
                Stop = obj.last_swap; 
            end

            if ~exist('RollingAvg','var')
                RollingAvg = 1; 
            end 
            energy_chain_to_plot = obj.energy_chain(1,:);
            plot_energy_chain(obj,energy_chain_to_plot,Start,RollingAvg,Stop)
       end 
       
       %Function to Simulate all fitted scenarios for the chosen samples.
       %Additionally, the data  will be stored in memory for easy plotting.
       function [obj] = simulate_fit_exp(obj)
           if isempty(obj.samples_chosen)
               obj.sample_patients;
               disp("No Sampled Parameters Found!") 
               disp(sprintf('Sampling %d times and assuming a converged start of %d iterations',obj.num_samples,obj.converged_start))
           elseif length(obj.samples_chosen) ~= obj.num_samples
                disp("Number of samples chosen does not equal to the number of samples set") 
                disp("Resampling to match sampled patients with number") 
                obj.sample_patients;
           elseif obj.current_converged_start ~= obj.converged_start
               disp("The minimum converged value changed, need to resample") 
               obj.sample_patients;
           else 
               disp("Reusing Previous Sampled Parameters!") 
               disp(sprintf('Sampled %d times and assumed a converged start of %d iterations',obj.num_samples,obj.converged_start))
           end 
           
           
           obj = prepare_empty_data_structure(obj); 

           for ith_sample = 1:obj.num_samples
               current_pars = obj.params_chain(1,:,obj.samples_chosen(ith_sample));
               obj.fitted_simulated_data.obv(:,:,:,ith_sample) = obj.key_struct.simulate_fit_trajectories(current_pars); 
               
               obj = calculate_and_store_psm(obj,ith_sample); 
           end
                
       end 
       
       %Function: Plot Single Analyte with all samples:
       function [] = plot_analyte(obj,Analyte,Expierment)
           %Check for valid inputs 
           if ~exist('Analyte','var')
               error("Need to provide analyte name")
           end 
           if ~exist('Expierment','var')
               Expierment_List = obj.key_struct.simulations_to_run
               error("Need to provide expierment name. Choose one of the above outputs from the Expierment_List shown above:") 
           end 
           
           obj = check_for_needed_resimulation(obj);
           
           
           data_to_plot = obj.get_simulated_data(Analyte,Expierment);
           
           %Plot the Data 
           figure 
           plot(obj.key_struct.simulation_time,data_to_plot)
           grid 
           xlabel("Time [Seconds]","FontSize",20)
           %set_xlabel(key_struct.time_information_plotting);
           set_ylabel(obj,Analyte,obj.key_struct.analyte_plotting_names);
           set_title(obj,Expierment,obj.key_struct.expierment_plotting_names);
               
       end 
       
       function [exported_simulated_data] = get_simulated_data(obj,Analyte,Expierment)
           %Check for valid inputs 
           if ~exist('Analyte','var')
               error("Need to provide analyte name")
           end 
           if ~exist('Expierment','var')
               Expierment_List = obj.key_struct.simulations_to_run
               error("Need to provide expierment name. Choose one of the above outputs from the Expierment_List shown above:") 
           end 
           
           %Determine the index of the analyte and expierment.First
           %objective is to determine if the analyte is a post simulation
           %modification or a direct observable
           index_analyte = find_index_of_analyte(obj,Analyte);
           
           %Find the index of the Expierment: 
           index_exp = find_index_of_expierment(obj,Expierment); 
           
           %Get the data: 
           if index_analyte > 0 
               exported_simulated_data = squeeze(obj.fitted_simulated_data.obv(:,index_analyte,index_exp,:));
           else 
               exported_simulated_data = squeeze(obj.fitted_simulated_data.post_modification.(Analyte).values(:,index_exp,:));
           end 
       end 
       
       function [] = plot_auto(obj)
           close all 
           obj = check_for_needed_resimulation(obj);
           
           for ith_fitted_analyte = 1:obj.key_struct.data.n_fitted_analytes
                analyte_i = obj.key_struct.data.analyte_names{ith_fitted_analyte};

                exp_data =  obj.key_struct.data.(analyte_i); 

                for ith_expierment = 1:length(exp_data.expierment_index) 
                    expierment_index = exp_data.expierment_index(ith_expierment);
                    expierment_name = obj.key_struct.simulations_to_run{expierment_index};
                    
                    obj.plot_analyte(analyte_i,expierment_name);
                    hold on 
                    errorbar(exp_data.time,exp_data.value(:,ith_expierment),... 
                             exp_data.sigma(:,ith_expierment),'LineStyle',"none",...
                             'Color','r','Marker','.','MarkerSize',20)%,...
                            %200,'r','filled')
                end
            end 
       end 
       
   end 
   methods(Hidden)
       %Intialization of the Plotting tool. 
       function obj = plotting_tool(LocOfKeyStruct,NumberOfSamplesWithOutReplacement,ConvergedStart)
           if nargin <1
               error("Must define the location of the ptemptest results")
           elseif nargin <2
               obj.loc_key_struct  = LocOfKeyStruct;
               obj.num_samples     = 100; 
               obj.converged_start = 1; 
           elseif nargin <3
               obj.loc_key_struct  = LocOfKeyStruct;
               obj.num_samples     = NumberOfSamplesWithOutReplacement; 
               obj.converged_start = 1;  
           else 
               obj.loc_key_struct = LocOfKeyStruct;
               obj.num_samples    = NumberOfSamplesWithOutReplacement; 
               obj.converged_start = ConvergedStart;
           end 
           
           obj = load_key_information(obj);
           
           %Set up fitted simulation data struct object
           obj.fitted_simulated_data.obv = [];
           obj.fitted_simulated_data.post_modification = [];
           
           %Set paths for simulation: 
           obj.key_struct.set_paths();
       end 
       
       %Loading the key_struct, parameter sets for each chain and swap, the
       %energy values for each chain and swap, and the last swap performed
       %in case the fit was finsihed early. 
       function obj = load_key_information(obj)
            temp_X = load(obj.loc_key_struct + "/key_struct.mat");

            progress_file_name = dir(temp_X.key_struct.output_location+"/"+temp_X.key_struct.progress_regex);
            additional = load(temp_X.key_struct.output_location+progress_file_name.name);

            obj.key_struct   = temp_X.key_struct; 
            obj.params_chain = additional.params_chain;
            obj.energy_chain = additional.energy_chain; 
            obj.last_swap    = additional.last_swap; 
       end 
       
       %Function to prepare the empty data structure for fitting
       function [obj] = prepare_empty_data_structure(obj)
           
           n_time   = length(obj.key_struct.simulation_time); 
           n_exp  = length(obj.key_struct.simulations_to_run);
           n_samp =  obj.num_samples;
           
           %Preparing the object that will store the observables directly
           n_obs    = obj.key_struct.n_obs;
           obj.fitted_simulated_data.obv =zeros(n_time,n_obs,n_exp,n_samp);
           
           %Preparing the object that will store the post simulation
           %modifications directly 
           temp_struct = struct();
           psm_names = fieldnames(obj.key_struct.post_simulation_analytes); 
           for ith_psm = 1:length(psm_names)
               psm_i = psm_names{ith_psm};
               temp_struct.(psm_i).values = zeros(n_time,n_exp,n_samp);
               temp_struct.(psm_i).parameters = zeros(n_samp,2);
           end 
           temp_struct.name_list = psm_names;
           obj.fitted_simulated_data.post_modification = temp_struct;                                 
       end 
       
       function [obj] = calculate_and_store_psm(obj,SampleN)
           
           %Calculate PSM 
           post_sim_modification = obj.key_struct.get_post_simulation_modifications(obj.fitted_simulated_data.obv(:,:,:,SampleN));
           
           %Store PSM 
           for ith_psm = 1:length(obj.fitted_simulated_data.post_modification.name_list)
               psm_i = obj.fitted_simulated_data.post_modification.name_list{ith_psm};
               obj.fitted_simulated_data.post_modification.(psm_i).values(:,:,SampleN)   = post_sim_modification.(psm_i).values;
               obj.fitted_simulated_data.post_modification.(psm_i).parameters(SampleN,:) = post_sim_modification.(psm_i).parameters;
           end 
           
       end 
       
       %Function to Plot the Energy Chains 
       function [] = plot_energy_chain(~,EnergyChain,Start,RollingAvg,Stop)
            figure 
            Y = movmean(EnergyChain,RollingAvg,2);
            Y = Y(:,Start:Stop);
            
            [num_chains,~] = size(Y);
            X = repmat(Start:Stop,num_chains,1)';
            
            plot(X,Y',"LineWidth",5)
            xlabel("Swap Iteration","FontSize",20)
            ylabel("Energy","FontSize",20)
            grid 
       end 
        
       
       %Find the index of the analyte. If the analyte is a direct
       %observable from BNGL then an index will be provided. If it is a
       %post simulation modification then -1 will be produced. 
       function [index_analyte] = find_index_of_analyte(obj,Analyte)
           if isfield(obj.key_struct.obsv_map,Analyte)
               %Observable Input:
               index_analyte = obj.key_struct.obsv_map.(Analyte);
           elseif isfield(obj.fitted_simulated_data.post_modification,Analyte)
               index_analyte = -1 ;
           else 
               %Future, produce a list of possible analyte
               error("No Analytes found, check name")
           end
           
       end 
       
       %Find the index of the expierment: 
       function [index_exp] = find_index_of_expierment(obj,Expierment)
           index_exp = strcmp(Expierment,obj.key_struct.simulations_to_run);
           
           if sum(index_exp)~=1 
               disp(obj.key_struct.simulations_to_run)
               error("Expierment Index Not Found. Check name")
           end 
           
       end 
       
       %Function that will check to see if re-simulation and re-sampling
       %needs to occur due to changes in the user defined variables. 
       function [obj] = check_for_needed_resimulation(obj)
           %Checking for existing simulated data. If empty, the simulation
           %will be performed. I will need to add an additional check to
           %ensure if re-sampling is done that the simulations are re-run. 
           if isempty(obj.fitted_simulated_data.obv)
               disp("Need to simulate and store data")
               obj.simulate_fit_exp
           elseif length(obj.samples_chosen) ~= obj.num_samples
               disp("Need to re-sample and simulate due to change in number of samples")
               obj.simulate_fit_exp
           elseif obj.current_converged_start ~= obj.converged_start
               disp("The minimum converged value changed, need to resimulate and resample") 
               obj.simulate_fit_exp
           end 
       end 
       
       %Functions to help set customizable figure labels: 
       function [NameOut] = check_for_user_defined_name(obj,Input,DictInput)
            if isfield(DictInput, Input)
                NameOut = DictInput.(Input);
            else 
                NameOut = Input ;
            end 
        end 

        function [] = set_ylabel(obj,Orginal,yLabelMap) 
            analyte_name = check_for_user_defined_name(obj,Orginal,yLabelMap);
            ylabel(analyte_name,"FontSize",20,'Interpreter','none')
        end 

        function [] = set_title(obj,Orginal,yLabelMap) 
            exp_name = check_for_user_defined_name(obj,Orginal,yLabelMap);
            title(exp_name,"FontSize",20,"Interpreter","none")
        end 

   end
end
