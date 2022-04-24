%% Impact of varying number of tasks on the objective function and number of replicas
%Preparing the data

N_min = 50;
N_max = 100;
N_stepSize = 5;

M_min = 5;
M_max = 20;
M_stepSize = 1;

number_of_simulations = 1;
checkConstraints = true;
dataObj = struct();
dataObj.rel_epsilon = 0.2;
guid = string(java.util.UUID.randomUUID.toString);
signature = strcat('[simulation_results1][N_min =  ', int2str(N_min), ', ', ...
    ' N_max =  ', int2str(N_max), ', ', ...
    ' N_stepSize =  ', int2str(N_stepSize), ', ', ...
    ' M_min =  ', int2str(M_min), ', ', ...
    ' M_max =  ', int2str(M_max),  ', ', ...
    ' M_stepSize =  ', int2str(M_stepSize),  ', ', ...
    ' n_sims =  ', int2str(number_of_simulations),  ']');
guid = strcat(signature, '_', guid);
save_to_file = true;

result = TREED_battery_aware_simulation(N_min, N_max, N_stepSize, M_min, M_max, M_stepSize, number_of_simulations, dataObj, checkConstraints);

N = ceil((N_max - N_min) ./ N_stepSize);
M = ceil((M_max - M_min) ./ M_stepSize);
n_vector = N_min:N_stepSize:N_max;
m_vector = M_min:M_stepSize:M_max
sim1 = struct();
for n=1:N
    for m = 1:M
        x = result{n,m}.all_sims{1}.x;
        temp_x = reshape(x, [m_vector(m), n_vector(n)]);
        num_of_replicas_per_task = sum(temp_x, 2);
        sim1.total_replicas(n, m) = sum(num_of_replicas_per_task);
        optimal_val = result{n,m}.all_sims{1}.optimalVal;
        sim1.all_recruited_workers(n, m) = sum(x);
        max_val = sum(result{n,m}.all_sims{1}.dataObj.objectiveFunction);
        max_val = max_val / n_vector(n);
        sim1.all_optimal_vals(n,m) = abs(optimal_val / max_val);
    end
end

if (save_to_file)
    save(strcat(guid, '.mat'));
end

%% Fixing number of workers and varying number of tasks
disp('Fixing number of workers and varying number of tasks');
i = size(sim1.all_optimal_vals, 1);
strresult = '';
tasks = 5;
for j = 1:size(sim1.all_optimal_vals, 2)
    opt_val = sim1.all_optimal_vals(i, j);
    
   strresult = strcat(strresult, '(', num2str(tasks),', ', num2str(opt_val), ')');
   tasks = tasks + M_stepSize;
end
disp(strresult);

%% Fixing number of tasks and varying number of workers
disp('Fixing number of tasks and varying number of workers');
j = size(sim1.all_optimal_vals, 2);
strresult = '';
workers = 50;
for i = 1:size(sim1.all_optimal_vals, 1)
   opt_val = sim1.all_optimal_vals(i, j);
   strresult = strcat(strresult, '(', num2str(workers),', ', num2str(opt_val), ')');
   workers = workers + N_stepSize;
end
disp(strresult);

%% Varying both the number of tasks and workers
disp('Varying both the number of tasks and workers');
strresult = '';
for i = 1:size(sim1.all_optimal_vals, 1)
    for j = 1:size(sim1.all_optimal_vals, 2)
        opt_val = sim1.all_optimal_vals(i, j);
        strresult = strcat(strresult, '(', num2str(tasks),', ', num2str(opt_val), ')');
    end
end
disp(strresult);
