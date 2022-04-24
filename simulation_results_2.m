%% Impact of varying number of tasks on the objective function and number of replicas
%Preparing the data
N_min = 100;
N_max = 100;
N_stepSize = 1;

M_min = 15;
M_max = 15;
M_stepSize = 1;

number_of_simulations = 1;
checkConstraints = true;

epsilons = 0.1:0.1:2;
 
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

dataObj.worker_hazzard_rate_fromval = 0.01;
dataObj.worker_hazzard_rate_toval = 0.5;
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
%dataObj.tasks_deadlines = round(dataObj.tasks_deadlines, 2);



num_of_replicas_per_task = cell(1, length(epsilons));
for eps = 1:length(epsilons)

    dataObj.rel_epsilon = epsilons(eps);
    result = TREED_battery_aware_simulation(N_min, N_max, N_stepSize, M_min, M_max, M_stepSize, number_of_simulations, dataObj, checkConstraints);
    x = result{1,1}.all_sims{1}.x;
    temp_x = reshape(x, [M_max, N_max]);
    num_of_replicas_per_task{eps} = sum(temp_x, 2);
    total_replicas(eps) = sum(num_of_replicas_per_task{eps});
    optimal_val = result{1,1}.all_sims{1}.optimalVal;
    all_recruited_workers(eps) = sum(x);
    max_val = sum(result{1,1}.all_sims{1}.dataObj.objectiveFunction);
    max_val = max_val / N_max;
    all_optimal_vals(eps) = abs(optimal_val / max_val);
end

strresult = '';
for i = 1:length(epsilons)
    strresult = strcat(strresult, '(', num2str(epsilons(i)),', ', num2str(total_replicas(i)), ')');
end
display(strresult);
