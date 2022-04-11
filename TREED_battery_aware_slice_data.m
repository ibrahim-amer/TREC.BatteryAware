function [dataObj] = TREED_battery_aware_slice_data(N, M, dataObj)
    dataObj.N = N;
    dataObj.M = M;
    dataObj.numOfVars = dataObj.N * dataObj.M + dataObj.M;
    dataObj.wrokers_freqs = dataObj.wrokers_freqs(:, 1:N);
    dataObj.workers_max_tasks = dataObj.workers_max_tasks(:, 1:N);
    dataObj.tasks_pdensity = dataObj.tasks_pdensity(:, 1:M);
    dataObj.tasks_dataSize = dataObj.tasks_dataSize(:, 1:M);
    dataObj.tasks_comp_delays = dataObj.tasks_comp_delays(:, 1:N * M);
    dataObj.tasks_deadlines = dataObj.tasks_deadlines(:, 1:M);
    dataObj.max_energies = dataObj.max_energies(:, 1:N);
    dataObj.comp_energies = dataObj.comp_energies(1:N*M, :);
    %            a     b    c   d   e   f_ub  f_lb  g_ub   g_lb
    dataObj.con_3b_size = N*M;
    dataObj.con_3c_size = N*M;
    dataObj.con_3d_size = N;
    dataObj.con_3e_size = M;
    dataObj.con_3f_size = 2;
    dataObj.con_3g_size = (N*M) + (N*M);
    dataObj.con_3h_size = M + M;
    A_b_n_rows= dataObj.con_3b_size + dataObj.con_3c_size + dataObj.con_3d_size + dataObj.con_3e_size + dataObj.con_3f_size + dataObj.con_3g_size + dataObj.con_3h_size;
    dataObj.A = dataObj.A(1:A_b_n_rows, :);
    dataObj.b = dataObj.b(1:A_b_n_rows, :);
    dataObj.Aeq = dataObj.Aeq(1:M, :);
    dataObj.beq = dataObj.beq(1:M, :);
    dataObj.ub = dataObj.ub(1, N * M + M);
    dataObj.lb = dataObj.lb(1, N * M + M);
    dataObj.operators = dataObj.operators(:, 1:A_b_n_rows);
end

