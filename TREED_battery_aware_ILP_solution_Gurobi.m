function result = TREED_battery_aware_ILP_solution_Gurobi(dataObj, checkConstraints, allow_debug)
    fprintf('#####[TREED_MILP_Solution_Gurobi] started!. N = %d, M = %d#####\n', dataObj.N, dataObj.M);
    %% Problem solution
    decVars = cell(1, dataObj.N * dataObj.M + dataObj.M);
    c1 = 1;
    %Prepare binary decision variables ,i.e., lambdas
    for i = 1:dataObj.N
        for j = 1:dataObj.M
            decVars{c1} = strcat('x', num2str(i), '_', num2str(j));
            c1 = c1 + 1;
        end
    end
    
    
    %decVars = {decvarBinary{:}, decvarCont{:}};
    %names = cat(1, decVars{:});
    names = reshape(decVars', [1, dataObj.numOfVars]);
    model.varnames = names;
    
    % Set objective:
    dataObj.objectiveFunction = [zeros(1, dataObj.N .* dataObj.M)];
    model.obj = dataObj.objectiveFunction;
    model.modelsense = 'min';
    
    % Linear constraints
    model.A = sparse(dataObj.A);
    model.rhs = dataObj.b;
    
    model.sense = dataObj.operators;

    
    % Decision variables types
    model.vtype = strcat(repmat('B', 1, dataObj.N .* dataObj.M), repmat('C', 1, dataObj.M));
       

    
    gurobi_write(model, 'mip1.lp');
    
   

    % Extract relevant Gurobi parameters from (subset of) options
    params = struct();
    %https://or.stackexchange.com/a/5166/5671
    params.IntegralityFocus = 1;
    minIntFeasTolVal = 1e-9;
    maxIntFeasTolVal = 1e-1;
    params.IntFeasTol = maxIntFeasTolVal;
%     params.NodefileStart = 0.5;
%     params.PreSparsify = 0;
    output = gurobi(model, params);
    result = struct();
    result.gurobiOutput = output;
    % Resolve model if status is INF_OR_UNBD
    if strcmp(output.status,'INFEASIBLE') || strcmp(output.status,'INF_OR_UNBD')
        params.DualReductions = 0;
        warning('Infeasible or unbounded, resolve without dual reductions to determine...');
        result.debug = struct();
        result.debug.iis = gurobi_iis(model);
        
%         penalties.lb = 2 .* ones(dataObj.numOfVars,1);
%         penalties.ub = 2 .* ones(dataObj.numOfVars,1);
%         penalties.rhs = 2 .* ones(length(dataObj.mcomb),1);
%         feasrelaxresult = gurobi_feasrelax(model, 0, false, penalties)
        result.x = [];
        result.optimalVal = [];
        result.status = output.status;
        if allow_debug
            result.debug.x = [];
            result.runtime = 0;
            
            result.debug.x = sym(zeros(dataObj.numOfVars, 1));
            for i = 1:dataObj.numOfVars
                sym_var = sym(names{i});
                assume(sym_var, 'real');
                result.debug.x(i, :) = sym_var;
            end

            result.debug.allConstraints = sym(zeros(size(dataObj.A, 1), 1));

            ctr = 1;
            for i=1:size(dataObj.A, 1)
                if dataObj.operators(1, ctr) == '<'
                    result.debug.allConstraints(i, :) = dataObj.A(i, :) * result.debug.x <= dataObj.b(i, :);
                elseif dataObj.operators(1, ctr) == '>'
                    result.debug.allConstraints(i, :) = dataObj.A(i, :) * result.debug.x >= dataObj.b(i, :);
                elseif dataObj.operators(1, ctr) == '='
                    result.debug.allConstraints(i, :) = dataObj.A(i, :) * result.debug.x == dataObj.b(i, :);
                end
                ctr = ctr + 1;
            end
        end
    elseif isfield(output, 'x')
        
        for j=1:dataObj.numOfVars
            fprintf('%s %e\n', names{j}, output.x(j));
        end
        fprintf('Obj: %e\n', output.objval);
        
        result.x = round(abs(output.x));
        result.runtime = output.runtime;
        result.status = output.status;
        result.x = result.x;
        result.stats = struct();
        result.optimalVal = output.objval;
        result.allConstraints = zeros(size(dataObj.A, 1), 1);
        result.debug = [];
        if checkConstraints
            ctr = 1;
            for i=1:size(dataObj.A, 1)
                if dataObj.operators(1, ctr) == '<'
                    result.allConstraints(i, :) = round(dataObj.A(i, :) * result.x, 3) <= dataObj.b(i, :);
                elseif dataObj.operators(1, ctr) == '>'
                    result.allConstraints(i, :) = round(dataObj.A(i, :) * result.x, 3) >= dataObj.b(i, :);
                elseif dataObj.operators(1, ctr) == '='
                    result.allConstraints(i, :) = round(dataObj.A(i, :) * result.x, 3) == dataObj.b(i, :);
                end
                ctr = ctr + 1;
            end
            result.allConstraintsSatisfied = sum(result.allConstraints) == length(result.allConstraints);
        end
    else 
        disp('The model is probablly infeasible or unbounded!');
        result.x = [];
        result.optimalVal = [];
        result.allConstraints = [];    
        result.status = 'INFEASIBLE';
    end
    fprintf('#####[MEC_Recruitment_MIP_Solution_Gurobi] finished!. N = %d, M = %d#####\n', dataObj.N, dataObj.M);
end