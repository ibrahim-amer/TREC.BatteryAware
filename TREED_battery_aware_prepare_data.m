function [dataObj] = TREED_battery_aware_prepare_data(dataObj)
    dataObj.numOfVars = dataObj.N .* dataObj.M + dataObj.M;
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
    %Constraint (3b)
    con_3b = zeros(dataObj.N * dataObj.M, dataObj.numOfVars);
    for i = 1:dataObj.N * dataObj.M
        row = zeros(1, dataObj.N .* dataObj.M + dataObj.M);
        row(i) = 1 .* dataObj.tasks_comp_delays(i);
        con_3b(i, :) = row;
    end
    con_3b_bounds = zeros(dataObj.N * dataObj.M, 1);
    ctr = 1;
    for i = 1:dataObj.N
        for j = 1:dataObj.M
            con_3b_bounds(ctr, :) = dataObj.tasks_deadlines(j);
            ctr = ctr + 1;
        end
    end
    dataObj.A = [dataObj.A; con_3b];
    dataObj.b = [dataObj.b; con_3b_bounds];
    dataObj.operators = [dataObj.operators repmat('<', 1, size(con_3b, 1))];
    % Constraint (3c) 
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
%     con_3b = zeros(dataObj.N * dataObj.M, dataObj.numOfVars);
%     for i = 1:dataObj.N * dataObj.M
%         row = zeros(1, dataObj.numOfVars);
%         row(i) = dataObj.comp_energies(i);
%         con_3b(i, :) = row;
%     end
    con_3c = zeros(dataObj.N * dataObj.M, dataObj.numOfVars);
    ctr = 1;
    for i = 0:dataObj.M:dataObj.N * dataObj.M - 1
        row = zeros(1, dataObj.numOfVars);
        for j = 1:dataObj.M
            row(1, i+j) = 1;            
        end
        con_3c(ctr, :) = row;
        ctr = ctr + 1;
    end
    
    con_3c_bounds = zeros(dataObj.N * dataObj.M, 1);
    ctr = 1;
    for i = 1:dataObj.N
        for j = 1:dataObj.M
            con_3c_bounds(ctr, :) = dataObj.max_energies(i);
            ctr = ctr + 1;
        end
    end
    
    dataObj.A = [dataObj.A; con_3c];
    dataObj.b = [dataObj.b;  con_3c_bounds];
    dataObj.operators = [dataObj.operators repmat('<', 1, size(con_3c, 1))];
    % Constraint (3d) (N)
    con_3d = zeros(dataObj.N, dataObj.numOfVars);
    ctr = 1;
    for i = 0:dataObj.M:dataObj.N * dataObj.M - 1
        row = zeros(1, dataObj.numOfVars);
        for j = 1:dataObj.M
            row(1, i+j) = 1;            
        end
        con_3d(ctr, :) = row;
        ctr = ctr + 1;
    end
    
    dataObj.A = [dataObj.A; con_3d];
    dataObj.b = [dataObj.b; dataObj.workers_max_tasks'];
    dataObj.operators = [dataObj.operators repmat('<', 1, size(con_3d, 1))];
    %Constraint (3e)
    con_3e = zeros(dataObj.M, dataObj.numOfVars);
    dataObj.Aeq = [];
    dataObj.beq = [];
    for i = 1:dataObj.M
        row = zeros(1, dataObj.numOfVars);
        count = 0;
        for j = 1:dataObj.N
            row(1, i+count) = 1;            
            count = count + dataObj.M;
        end
        max_delay = zeros(1, dataObj.N .* dataObj.M + dataObj.M);
        max_delay(1, dataObj.N .* dataObj.M + i) = -1;
        row = row + max_delay;
        con_3e(i, :) = row;
    end
    dataObj.Aeq = [dataObj.Aeq; con_3e];
    dataObj.beq = [dataObj.beq; zeros(dataObj.M, 1);];
    %merge inequalities with equalities 
    dataObj.A = [dataObj.A; dataObj.Aeq];
    dataObj.b = [dataObj.b; dataObj.beq];
    dataObj.operators = [dataObj.operators repmat('=', 1, size(dataObj.Aeq, 1))];
    %Constraint (3f) RHS
    con_3f_rhs = zeros(1, dataObj.N .* dataObj.M + dataObj.M);
    con_3f_rhs(1, (dataObj.N * dataObj.M + 1) : end) = 1;
    dataObj.A = [dataObj.A; con_3f_rhs];
    dataObj.b = [dataObj.b; sum(dataObj.workers_max_tasks)];
    dataObj.operators = [dataObj.operators repmat('<', 1, size(con_3f_rhs, 1))];
    %Constraint (3f) LHS
    con_3f_lhs = zeros(1, dataObj.N .* dataObj.M + dataObj.M);
    con_3f_lhs(1, (dataObj.N * dataObj.M + 1) : end) = -1;
    dataObj.A = [dataObj.A; con_3f_lhs];
    dataObj.b = [dataObj.b; -1 .* (dataObj.M ./ dataObj.N)];
    dataObj.operators = [dataObj.operators repmat('<', 1, size(con_3f_lhs, 1))];
    
    dataObj.ub = zeros(1, dataObj.N * dataObj.M + dataObj.M);
    dataObj.lb = zeros(1, dataObj.N * dataObj.M + dataObj.M);
    dataObj.ub(1, 1:dataObj.N * dataObj.M) = 1;
    dataObj.lb(1, 1:dataObj.N * dataObj.M) = 0; %already satisfied
    dataObj.ub(1, dataObj.N * dataObj.M + 1:end) = dataObj.N;
    dataObj.lb(1, dataObj.N * dataObj.M + 1: end) = 1;
    
    
    %Upper and lower bounds
    %lambdas upper bounds
    lambdasBounds = diag([ones(dataObj.N .* dataObj.M, 1); zeros(dataObj.M, 1)]);
    lambdasBounds = lambdasBounds(1:dataObj.N * dataObj.M, :);
    dataObj.A = [dataObj.A; lambdasBounds];
    dataObj.b = [dataObj.b; dataObj.ub(1, 1:dataObj.N * dataObj.M)'];
    dataObj.operators = [dataObj.operators repmat('<', 1, dataObj.N * dataObj.M)];
    %lambdas lower bounds
    dataObj.A = [dataObj.A; lambdasBounds];
    dataObj.b = [dataObj.b; dataObj.lb(1, 1:dataObj.N * dataObj.M)'];
    dataObj.operators = [dataObj.operators repmat('>', 1, dataObj.N * dataObj.M)];
    %replicas upper bounds
    replicasBounds = diag([zeros(dataObj.N .* dataObj.M, 1); ones(dataObj.M, 1)]);
    replicasBounds = replicasBounds(dataObj.N * dataObj.M + 1:end, :);
    dataObj.A = [dataObj.A; replicasBounds];
    dataObj.b = [dataObj.b; dataObj.ub(1, dataObj.N * dataObj.M + 1:end)'];
    dataObj.operators = [dataObj.operators repmat('<', 1, dataObj.M)];
    %replicas lower bounds
    dataObj.A = [dataObj.A; replicasBounds];
    dataObj.b = [dataObj.b; dataObj.lb(1, dataObj.N * dataObj.M + 1:end)'];
    dataObj.operators = [dataObj.operators repmat('>', 1, dataObj.M)];
    
    %%
end

