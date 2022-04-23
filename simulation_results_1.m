%% Impact of varying number of tasks on the objective function and number of replicas
%Preparing the data
N_min = 50;
N_max = 100;
N_stepSize = 10;

M_min = 5;
M_max = 10;
M_stepSize = 1;

number_of_simulations = 1;
checkConstraints = true;
dataObj = struct();
dataObj.rel_epsilon = 0.2;
result = TREED_battery_aware_simulation(N_min, N_max, N_stepSize, M_min, M_max, M_stepSize, number_of_simulations, dataObj, checkConstraints);

N = ceil((N_max - N_min) ./ N_stepSize);
M = ceil((M_max - M_min) ./ M_stepSize);
n_vector = linspace(N_min, N_max, N_stepSize);
m_vector = linspace(M_min, M_max, M_stepSize);
for n=1:N
    for m = 1:M
        x = result{n,m}.all_sims{1}.x;
        temp_x = reshape(x, [N, M]);
        num_of_replicas_per_task = sum(temp_x, 1);
        optimal_val = result{n,m}.all_sims{1}.optimalVal;
        all_recruited_workers(n, m) = sum(x);
        all_optimal_vals(n,m) = optimal_val;
    end
end