function [result] = TREED_battery_aware_simulation(N_min, N_max, N_stepSize, M_min, M_max, M_stepSize, number_of_simulations, checkConstraints)
%%
arguments
    N_min (1, 1) double = 2
    N_max (1, 1) double = 100
    N_stepSize (1, 1) double = 2
    M_min (1, 1) double = 2
    M_max (1, 1) double = 20
    M_stepSize (1, 1) = 2
    number_of_simulations (1, 1) double = 1.0
    checkConstraints (1, 1) logical = 1
end
    save_to_file = false;
    %%Preamble
    guid = string(java.util.UUID.randomUUID.toString);
    signature = strcat('[scenario1][N_min =  ', int2str(N_min), ', ', ...
        ' N_max =  ', int2str(N_max), ', ', ...
        ' N_stepSize =  ', int2str(N_stepSize), ', ', ...
        ' M_min =  ', int2str(M_min), ', ', ...
        ' M_max =  ', int2str(M_max),  ', ', ...
        ' M_stepSize =  ', int2str(M_stepSize),  ', ', ...
        ' n_sims =  ', int2str(number_of_simulations),  ']');
    guid = strcat(signature, '_', guid);
    allResults_size = (ceil((N_max - N_min + 1) ./ N_stepSize) .* ceil((M_max - M_min + 1) ./ M_stepSize)) + 1;
    allResults = cell(1, allResults_size);
    ctr = 1;
    %%Preparing data
    dataObj.N = N_max;
    dataObj.M = M_max;
    dataObj = TREED_battery_aware_prepare_data(dataObj);
    max_iter_random_policy = 500;
    allow_debug = true;
    
    %%Run simulations
    for N = N_min:N_stepSize:N_max
        for M = M_min:M_stepSize:M_max
            fprintf('###########################################################\n');
            fprintf(strcat('N = ', int2str(N), ' ', ' OUT OF N = ', int2str(N_max), ' M = ', int2str(M), ...
                    ' ', ' OUT OF M = ', int2str(M_max), ' STARTS!!!', '\n'));
            result = struct();
            result.results = cell(1, number_of_simulations);
            result.stats = struct();
            %dataObj = TREED_battery_aware_slice_data(N, M, dataObj);
            
            %result.random_policy = TREED_maximize_replicas_random_policy(dataObj, max_iter_random_policy);
            %result.baseline = TREED_maximize_replicas_baseline(dataObj);
           
            result.stats.averageOptimalVal = double(0);
            result.stats.averageRuntime = double(0);
            for sim = 1:number_of_simulations
                result.results{sim} = TREED_battery_aware_ILP_solution_Gurobi(dataObj, checkConstraints, allow_debug);
                
                if (strcmp(result.results{sim}.status, 'OPTIMAL') || strcmp(result.results{sim}.status, 'INTERRUPTED'))
                    result.stats.averageOptimalVal = result.stats.averageOptimalVal + result.results{sim}.optimalVal;
                    result.stats.averageRuntime = result.stats.averageRuntime + result.results{sim}.runtime;
                end
            end
            if (strcmp(result.results{sim}.status, 'OPTIMAL') || strcmp(result.results{sim}.status, 'INTERRUPTED'))
                result.stats.averageOptimalVal = result.stats.averageOptimalVal ./ number_of_simulations;
                result.stats.averageRuntime = result.stats.averageRuntime ./ number_of_simulations;
            end
            
            
            
            allResults{ctr} = result;
            ctr = ctr + 1;
            if (save_to_file)
                save(strcat(guid, '.mat'));
            end
            result.dataObj = dataObj;
            %clc;
            fprintf(strcat('N = ', int2str(N), ' ', ' OUT OF N = ', int2str(N_max), ...
                ' M = ', int2str(M), ' ', ' OUT OF K = ', int2str(M_max), ' ENDS!!!', '\n'));
            fprintf('###########################################################\n');
        end
    end
    if (save_to_file)
        save(strcat(guid, '.mat'));
    end
    
end

