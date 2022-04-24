%% Comparing TRUE approach with a method that only maximizes the number of replicas
%Preparing the data
N_min = 50;
N_max = 100;
N_stepSize = 10;

M_min = 5;
M_max = 30;
M_stepSize = 1;

number_of_simulations = 1;
checkConstraints = true;
dataObj.max_energy = 2.8;%was 1
 
dataObj = struct();
dataObj.N = N_max;
dataObj.M = M_max;
dataObj.numOfVars = dataObj.N * dataObj.M;
%Communication model parameters
    
dataObj.trans_power = 50 * 1e-3; %50 mWatt;
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

dataObj.rel_epsilon = 0.8;
result1 = TREED_battery_aware_simulation(N_min, N_max, N_stepSize, M_min, M_max, M_stepSize, number_of_simulations, dataObj, checkConstraints);
result2 = TREED_maximize_replicas_simulation(N_min, N_max, N_stepSize, M_min, M_max, M_stepSize, number_of_simulations, dataObj, checkConstraints);
N = ceil((N_max - N_min + 1) ./ N_stepSize);
M = ceil((M_max - M_min + 1) ./ M_stepSize);
n_vector = N_min:N_stepSize:N_max;
m_vector = M_min:M_stepSize:M_max;
sim4 = struct();
sim4.result1 = struct();
sim4.result2 = struct();
strresult1 = '';
strresult2 = '';
for n=1:N
    for m = 1:M
        try
            x = result1{n,m}.all_sims{1}.x;
            sim4.result1.rel_scores(n, m) = result1{n,m}.dataObj.workers_tasks_rel_prop * x;
            sim4.result1.Ms(n,m) = result1{n,m}.dataObj.M;
            x_axis = sim4.result1.Ms(n,m);
            y_axis = sim4.result1.rel_scores(n, m);
            
            strresult1 = strcat(strresult1, '(', num2str(x_axis),', ', num2str(y_axis), ')');

            x = result2{n,m}.all_sims{1}.x;
            sim4.result2.Ms(n,m) = result2{n,m}.dataObj.M;
            sim4.result2.Ns(n,m) = result2{n,m}.dataObj.N;
            sim4.result2.rel_scores(n, m) = result2{n,m}.dataObj.workers_tasks_rel_prop * x(1:sim4.result2.Ms(n,m) * sim4.result2.Ns(n,m), :);
            
            
            x_axis = sim4.result2.Ms(n,m);
            y_axis = sim4.result2.rel_scores(n, m);
            
            strresult2 = strcat(strresult2, '(', num2str(x_axis),', ',  num2str(y_axis), ')');
        catch me
            y = 1;
        end
    end
end
%% Varying number of tasks and displaying the results against the reliability score
strresult1 = '';
for n = N:N
    for j = 1:M
        x_axis = sim4.result1.Ms(n,m);
        y_axis = sim4.result1.rel_scores(n, m);
        strresult1 = strcat(strresult1, '(', num2str(x_axis),', ', num2str(y_axis), ')');
    end
end
disp('%% Varying number of tasks and displaying the results against the reliability score');
disp(strresult1);
%% Varying number of workers and displaying the results against the reliability score
strresult1 = '';
for n = 1:N
    for j = M:M
        x_axis = sim4.result2.Ms(n,m);
        y_axis = sim4.result2.rel_scores(n, m);
        strresult1 = strcat(strresult1, '(', num2str(x_axis),', ', num2str(y_axis), ')');
    end
end
disp('%% Varying number of workers and displaying the results against the reliability score');
disp(strresult1);

