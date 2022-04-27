%% Impact of varying number of tasks on on the number of completed tasks, between TRUE, TRUE-RD, and MTREED
%Preparing the data
N_min = 100;
N_max = 100;
N_stepSize = 1;

M_min = 1;
M_max = 10;
M_stepSize = 1;

number_of_simulations = 1;
checkConstraints = true;


 
dataObj = struct();
dataObj.N = N_max;
dataObj.M = M_max;
dataObj.numOfVars = dataObj.N * dataObj.M;
%Communication model parameters
    
dataObj.trans_power = 50 * 1e-3 %50 mWatt;
dataObj.path_loss_exp = 2;
dataObj.sigma_sq = 1e-11;
dataObj.controller_bandwidth = 10e6;

%Reliablity level threshold


dataObj.worker_CPU_FromVal = 1e9;
dataObj.worker_CPU_ToVal = 4e9;
dataObj.workers_freqs = dataObj.worker_CPU_FromVal + (dataObj.worker_CPU_ToVal - dataObj.worker_CPU_FromVal) * rand(1, dataObj.N);  % size  = N
%dataObj.workers_freqs = round(dataObj.workers_freqs, 2);

%Workers maximum number of tasks

dataObj.worker_max_tasks_fromval = 1;
dataObj.worker_max_tasks_toval = 3;
dataObj.workers_max_tasks = dataObj.worker_max_tasks_fromval + (dataObj.worker_max_tasks_toval - dataObj.worker_max_tasks_fromval) * rand(1, dataObj.N);  % size  = N
dataObj.workers_max_tasks = round(dataObj.workers_max_tasks);

%Workers distance from the controller

dataObj.worker_distances_fromval = 5;
dataObj.worker_distances_toval = 50;
dataObj.workers_distances = dataObj.worker_distances_fromval + (dataObj.worker_distances_toval - dataObj.worker_distances_fromval) * rand(1, dataObj.N);  % size  = N
%dataObj.workers_distances = round(dataObj.workers_distances);

%Workers Rayleigh coefficient

dataObj.workers_rayleigh = exprnd(1, [1, dataObj.N]); %mu = 1 -->unit mean

%Workers hazard rates

dataObj.worker_hazzard_rate_fromval = 0.08;
dataObj.worker_hazzard_rate_toval = 0.7;
dataObj.workers_hazard_rates = dataObj.worker_hazzard_rate_fromval + (dataObj.worker_hazzard_rate_toval - dataObj.worker_hazzard_rate_fromval) * rand(1, dataObj.N);  % size  = N
%dataObj.workers_hazard_rates = round(dataObj.workers_hazard_rates);

%Tasks' Processing Density

dataObj.task_pdensity_fromVal = 1e2;
dataObj.task_pdensity_toVal = 5e2;
dataObj.tasks_pdensity = dataObj.task_pdensity_fromVal + (dataObj.task_pdensity_toVal - dataObj.task_pdensity_fromVal) * rand(1, dataObj.M);  % size  = M
%dataObj.tasks_pdensity = round(dataObj.tasks_pdensity, 2);


%Tasks data size

dataObj.task_dataSize_fromVal = 1e6;
dataObj.task_dataSize_toVal = 20e6;
dataObj.tasks_dataSize = dataObj.task_dataSize_fromVal + (dataObj.task_dataSize_toVal - dataObj.task_dataSize_fromVal) * rand(1, dataObj.M);  % size  = M
%dataObj.tasks_dataSize = round(dataObj.tasks_dataSize, 2);


%Tasks CPU requirement

dataObj.task_CPU_fromVal = 1e9;
dataObj.task_CPU_toVal = 2e9;
dataObj.tasks_CPU_req = dataObj.task_CPU_fromVal + (dataObj.task_CPU_toVal - dataObj.task_CPU_fromVal) * rand(1, dataObj.M);  % size  = M
%dataObj.tasks_CPU_req = round(dataObj.tasks_dataSize, 2);


%Tasks deadlines - uniformly distributed

dataObj.task_deadline_fromVal = 8;%was 5
dataObj.task_deadline_toVal = 10;%was 20
dataObj.tasks_deadlines = dataObj.task_deadline_fromVal + (dataObj.task_deadline_toVal - dataObj.task_deadline_fromVal) * rand(1, dataObj.M);  % size  = M
dataObj.delay_dividend = 2;
%dataObj.tasks_deadlines = round(dataObj.tasks_deadlines, 2);

N = ceil((N_max - N_min + 1) ./ N_stepSize);
M = ceil((M_max - M_min + 1) ./ M_stepSize);
n_vector = N_min:N_stepSize:N_max;
m_vector = M_min:M_stepSize:M_max;

epsilons =  0.8:0.2:1.2;
all_results = cell(1, length(epsilons));
for eps = 1:length(epsilons)

    dataObj.rel_epsilon = epsilons(eps);
    
    result1 = TREED_battery_aware_simulation(N_min, N_max, N_stepSize, M_min, M_max, M_stepSize, number_of_simulations, dataObj, checkConstraints);
    result2 = TREED_battery_aware_simulation_rel_dropped(N_min, N_max, N_stepSize, M_min, M_max, M_stepSize, number_of_simulations, dataObj, checkConstraints);
    result3 = TREED_maximize_replicas_simulation(N_min, N_max, N_stepSize, M_min, M_max, M_stepSize, number_of_simulations, dataObj, checkConstraints);
    
    all_results{eps}.result1 = result1;
    all_results{eps}.result2 = result2;
    all_results{eps}.result3 = result3;
    
end

plot1 = '';
plot2 = '';
plot3 = '';
for eps = 1:length(epsilons)
    for m = 1:M
        try
            failure_percentage = 0.5;
            x1 = all_results{eps}.result1{1,m}.all_sims{1}.x;
            x2 = all_results{eps}.result2{1,m}.all_sims{1}.x;
            x3 = all_results{eps}.result3{1,m}.all_sims{1}.x;
            x3 = x3(1:(length(x3) - current_m));
            current_m = all_results{eps}.result1{1,m}.dataObj.M;
            
            x1_reshaped = reshape(x1, [current_m, N_max]);
            x2_reshaped = reshape(x2, [current_m, N_max]);
            x3_reshaped = reshape(x3, [current_m, N_max]);
            
            true_replicas = sum(x1_reshaped, 2);
            true_RD_replicas = sum(x2_reshaped, 2);
            mtreed_replicas = sum(x3_reshaped, 2);
            
            
            
            % Calculate failure probability
            failure_probs = 1 - all_results{eps}.result1{1,m}.dataObj.workers_tasks_rel_prop;
            true_tasks_status = (x1_reshaped .* failure_probs) > failure_percentage;
            true_RD_tasks_status = (x2_reshaped .* failure_probs) > failure_percentage;
            mtreed_tasks_status = (x3_reshaped .* failure_probs) > failure_percentage;
            
            true_nodes_failed = sum(true_tasks_status, 2);
            true_RD_nodes_failed = sum(true_RD_tasks_status, 2);
            mtreed_nodes_failed = sum(mtreed_tasks_status, 2);
            
            
            
            true_percentage = sum((failure_probs' .* x1) > failure_percentage) ./ current_m;
            true_RD_percentage = sum((failure_probs' .* x2) > failure_percentage) ./ current_m;
            mtreed_percentage = sum((failure_probs' .* x3) > failure_percentage) ./ current_m;

            plot1 = strcat(plot1, '(', num2str(current_m),', ', num2str(true_percentage), ')');
            plot2 = strcat(plot2, '(', num2str(current_m),', ', num2str(true_RD_percentage), ')');
            plot3 = strcat(plot3, '(', num2str(current_m),', ', num2str(mtreed_percentage), ')');
        catch me
            disp('error occured!');
        end
    end
disp(strcat('Percentage vs #Tasks, for epsilon = ', num2str(epsilons(eps))));
disp(plot1);
disp(plot2);
disp(plot3);
disp('#########################################################');
end


