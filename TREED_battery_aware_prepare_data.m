function [dataObj] = TREED_battery_aware_prepare_data(dataObj)
    dataObj.numOfVars = dataObj.N .* dataObj.M;
    %% Problem formulation and system model can be found here: 
    % https://skillful-honesty-f66.notion.site/Meeting-Preparation-March-23-2022-b7f0da29e5694554ba0f07d0acefe679
    %% Preparing Constants: 
    % density: 100MHz - 300 MHz
    % CPU frequency: 1-4 GhZ
    % data size: 1-50 MB
    %%Preparing constants
    %Kappa = 1e-29
    %Worker CPU Frequency
    dataObj.worker_CPU_FromVal = 1e9;
    dataObj.worker_CPU_ToVal = 4e9;
    dataObj.wrokers_freqs = dataObj.worker_CPU_FromVal + (dataObj.worker_CPU_ToVal - dataObj.worker_CPU_FromVal) * rand(1, dataObj.N);  % size  = N
    dataObj.wrokers_freqs = round(dataObj.wrokers_freqs, 2);
    
    %Workers maximum number of tasks
    dataObj.worker_max_tasks_fromval = 1;
    dataObj.worker_max_tasks_toval = 3;
    dataObj.workers_max_tasks = dataObj.worker_max_tasks_fromval + (dataObj.worker_max_tasks_toval - dataObj.worker_max_tasks_fromval) * rand(1, dataObj.N);  % size  = N
    dataObj.workers_max_tasks = round(dataObj.workers_max_tasks);
    
    %Workers battery capacities
    %Workers maximum number of tasks
    dataObj.worker_battery_fromval = 0;
    dataObj.worker_battery_toval = 1;
    dataObj.workers_battery_caps = dataObj.worker_battery_fromval + (dataObj.worker_battery_toval - dataObj.worker_battery_fromval) * rand(1, dataObj.N);  % size  = N
    %dataObj.workers_battery_caps = round(dataObj.workers_max_tasks);
    dataObj.workers_max_batt_cap = max(dataObj.workers_battery_cap);
    %Kappa 
    dataObj.kappa = 1e-29;
    
    %Tasks' Processing Density
    dataObj.task_pdensity_fromVal = 1e2;
    dataObj.task_pdensity_toVal = 10e2;
    dataObj.tasks_pdensity = dataObj.task_pdensity_fromVal + (dataObj.task_pdensity_toVal - dataObj.task_pdensity_fromVal) * rand(1, dataObj.M);  % size  = M
    dataObj.tasks_pdensity = round(dataObj.tasks_pdensity, 2);
    
    %Tasks data size
    dataObj.task_dataSize_fromVal = 1e6;
    dataObj.task_dataSize_toVal = 50e6;
    dataObj.tasks_dataSize = dataObj.task_dataSize_fromVal + (dataObj.task_dataSize_toVal - dataObj.task_dataSize_fromVal) * rand(1, dataObj.M);  % size  = M
    dataObj.tasks_dataSize = round(dataObj.tasks_dataSize, 2);
    
    %Tasks CPU requirement
    dataObj.task_CPU_fromVal = 1;
    dataObj.task_CPU_toVal = 2;
    dataObj.tasks_CPU_req = dataObj.task_CPU_fromVal + (dataObj.task_CPU_toVal - dataObj.task_CPU_fromVal) * rand(1, dataObj.M);  % size  = M
    dataObj.tasks_CPU_req = round(dataObj.tasks_dataSize, 2);
    
    %Computation delays
    tasks_specs = dataObj.tasks_pdensity .* dataObj.tasks_dataSize; % vectorSize = M
    dataObj.tasks_comp_delays = [];
    for i=1:dataObj.N
        dataObj.tasks_comp_delays = [dataObj.tasks_comp_delays (tasks_specs ./ dataObj.wrokers_freqs(i))];
    end
    %max_delay = task_dataSize_toVal * task_pdensity_toVal / worker_CPU_FromVal;
    %dataObj.tasks_comp_delays = dataObj.tasks_comp_delays / max_delay;
    
    %Tasks deadlines - uniformly distributed
    dataObj.task_deadline_fromVal = 20;%was 5
    dataObj.task_deadline_toVal = 30;%was 20
    dataObj.tasks_deadlines = dataObj.task_deadline_fromVal + (dataObj.task_deadline_toVal - dataObj.task_deadline_fromVal) * rand(1, dataObj.M);  % size  = M
    dataObj.tasks_deadlines = round(dataObj.tasks_deadlines, 2);
    
    %% Preparing constraints matrix A 
    dataObj.A = [];
    dataObj.b = [];
    dataObj.operators = [];
    %% Constraint w^{\text{CPU}}_i\lambda_{ij} >= t^{\text{CPU}}_j (b)
    con_b = zeros(dataObj.N * dataObj.M, dataObj.numOfVars);
    con_b_bounds = zeros(dataObj.N * dataObj.M, 1);
    ctr = 1;
    for i = 1:dataObj.N
        for j = 1:dataObj.M
            row = zeros(1, dataObj.numOfVars);
            row(i + j - 1) = 1 .* dataObj.wrokers_freqs(i);
            con_b(ctr, :) = row;
            con_b_bounds(ctr, :) = dataObj.tasks_CPU_req(j);
        end
    end
    dataObj.A = [dataObj.A; con_b];
    dataObj.b = [dataObj.b; con_b_bounds];
    dataObj.operators = [dataObj.operators repmat('>', 1, size(con_b, 1))];
    dataObj.con_b_size = size(con_b, 1);
    
    
    %% Constraint \sum^N_{i =1} w^{\text{cap}}_i\lambda_{ij} >= w^{\text{cap}}_{\max} (c)
    con_c = zeros(dataObj.M, dataObj.numOfVars);
    con_c_bounds = zeros(dataObj.M, 1);
    ctr = 1;
    for j = 1:dataObj.M
        row = zeros(1, dataObj.numOfVars);
        for i = j:dataObj.M:dataObj.N * dataObj.M
            row(i) = dataObj.workers_battery_caps(i / dataObj.max_batt_cap);
            con_c_bounds(ctr) = dataObj.workers_max_batt_cap;
        end
        con_c(j) = row;
    end
    dataObj.A = [dataObj.A; con_c];
    dataObj.b = [dataObj.b; con_c_bounds];
    dataObj.operators = [dataObj.operators repmat('>', 1, size(con_c, 1))];
    dataObj.con_c_size = size(con_c, 1);
    
    %% Constraint \tau^{\text{comp}}_{ij}\lambda_{ij} <=  t^{\text{deadline}}_j (d)
    con_d = zeros(dataObj.N * dataObj.M, dataObj.numOfVars);
    for i = 1:dataObj.N * dataObj.M
        row = zeros(1, dataObj.N .* dataObj.M + dataObj.M);
        row(i) = 1 .* dataObj.tasks_comp_delays(i);
        con_d(i, :) = row;
    end
    con_d_bounds = zeros(dataObj.N * dataObj.M, 1);
    ctr = 1;
    for i = 1:dataObj.N
        for j = 1:dataObj.M
            con_d_bounds(ctr, :) = dataObj.tasks_deadlines(j);
            ctr = ctr + 1;
        end
    end
    dataObj.A = [dataObj.A; con_d];
    dataObj.b = [dataObj.b; con_d_bounds];
    dataObj.operators = [dataObj.operators repmat('<', 1, size(con_d, 1))];
    dataObj.con_d_size = size(con_d, 1);
    
    %% Constraint \sum^N_{i = 1}\lambda_{ij} >= 1 (e)
    con_e = zeros(dataObj.M, dataObj.numOfVars);
    for j = 1:dataObj.M
        row = zeros(1, dataObj.numOfVars);
        for ij = 0:dataObj.M:dataObj.N * dataObj.M - 1
            row(j + ij) = 1; 
        end
        con_e(j) = row;
    end
    dataObj.A = [dataObj.A; con_e];
    dataObj.b = [dataObj.b; ones(size(con_e, 1), 1)];
    dataObj.operators = [dataObj.operators repmat('>', 1, size(con_e, 1))];
    dataObj.con_e_size = size(con_e, 1);
    
    %% Constraint \sum^{M}_{j=1} \lambda_{ij} <= w^{\text{tasks}}_i(f)
    con_f = zeros(dataObj.N, dataObj.numOfVars);
    ctr = 1;
    for i = 0:dataObj.M:dataObj.N * dataObj.M - 1
        row = zeros(1, dataObj.numOfVars);
        for j = 1:dataObj.M
            row(1, i+j) = 1;            
        end
        con_f(ctr, :) = row;
        ctr = ctr + 1;
    end
    
    dataObj.A = [dataObj.A; con_f];
    dataObj.b = [dataObj.b; dataObj.workers_max_tasks'];
    dataObj.operators = [dataObj.operators repmat('<', 1, size(con_f, 1))];
    dataObj.con_f_size = size(con_f, 1);
    
    %% Constraint (ee) 
    dataObj.max_energy = 3;%was 1
    dataObj.max_energies = dataObj.max_energy .* ones(1, dataObj.N);%size N
    energy_task_specs_tmp = dataObj.kappa .* dataObj.tasks_dataSize .* dataObj.tasks_pdensity; %size M
    dataObj.comp_energies = zeros(dataObj.N * dataObj.M, 1); %size M * N, 1
    ctr = 1;
    for i = 1:dataObj.N
        for j = 1:dataObj.M
            dataObj.comp_energies(ctr) = energy_task_specs_tmp(j) .* dataObj.wrokers_freqs(i)^2;
            ctr = ctr + 1;
        end
    end
    con_ee = zeros(dataObj.N * dataObj.M, dataObj.numOfVars);
    ctr = 1;
    for i = 0:dataObj.M:dataObj.N * dataObj.M - 1
        row = zeros(1, dataObj.numOfVars);
        for j = 1:dataObj.M
            row(1, i+j) = 1;            
        end
        con_ee(ctr, :) = row;
        ctr = ctr + 1;
    end
    
    con_ee_bounds = zeros(dataObj.N * dataObj.M, 1);
    ctr = 1;
    for i = 1:dataObj.N
        for j = 1:dataObj.M
            con_ee_bounds(ctr, :) = dataObj.max_energies(i);
            ctr = ctr + 1;
        end
    end
    
    dataObj.A = [dataObj.A; con_ee];
    dataObj.b = [dataObj.b;  con_ee_bounds];
    dataObj.operators = [dataObj.operators repmat('<', 1, size(con_ee, 1))];
    dataObj.con_e_size = size(con_ee, 1);

    %% Constraint (g) RHS
    
    dataObj.ub = zeros(1, dataObj.N * dataObj.M);
    dataObj.lb = zeros(1, dataObj.N * dataObj.M);
    dataObj.ub(1, 1:dataObj.N * dataObj.M) = 1;
    dataObj.lb(1, 1:dataObj.N * dataObj.M) = 0; %already satisfied
    
    
    %Upper and lower bounds
    %lambdas upper bounds
    lambdasBounds = diag([ones(dataObj.N .* dataObj.M, 1);]);
    lambdasBounds = lambdasBounds(1:dataObj.N * dataObj.M, :);
    dataObj.A = [dataObj.A; lambdasBounds];
    dataObj.b = [dataObj.b; dataObj.ub(1, 1:dataObj.N * dataObj.M)'];
    dataObj.operators = [dataObj.operators repmat('<', 1, dataObj.N * dataObj.M)];
    %lambdas lower bounds
    dataObj.A = [dataObj.A; lambdasBounds];
    dataObj.b = [dataObj.b; dataObj.lb(1, 1:dataObj.N * dataObj.M)'];
    dataObj.operators = [dataObj.operators repmat('>', 1, dataObj.N * dataObj.M)];
    dataObj.bounds_size = size(2 * lambdasBounds, 1);
    
    %%
end

