function  [x_optimal, cash_optimal] = strat_min_variance(x_init, cash_init, mu, Q, cur_prices)
    addpath('/Applications/CPLEX_Studio128/cplex/matlab/x86-64_osx');
    
    % number of stocks
    n = 20;
    
    % Optimization problem data
    lb = zeros(n,1);
    ub = inf*ones(n,1);
    A  = ones(1,n);
    b  = 1;
    
    % Compute minimum variance portfolio
    cplex1 = Cplex('min_Variance');
    cplex1.addCols(zeros(n,1), [], lb, ub);
    cplex1.addRows(b, A, b);
    cplex1.Model.Q = 2*Q;
    cplex1.Param.qpmethod.Cur = 6; % concurrent algorithm
    cplex1.Param.barrier.crossover.Cur = 1; % enable crossover
    cplex1.DisplayFunc = []; % disable output to screen
    cplex1.solve();
    
    % Display minimum variance portfolio
    w_minVar = cplex1.Solution.x; % weight for each stock  20*1
    cur_money = cur_prices * x_init + cash_init; % scalar
    allocated_money = w_minVar * cur_money; % 20*1
    x_optimal = round(allocated_money' ./ cur_prices); % 1*20
    transaction = abs((x_init' - x_optimal)) * cur_prices' * 0.005; %scalar
    cash_optimal = cur_money - cur_prices * x_optimal' - transaction; %scalar
    
    if cash_optimal < 0
        % Find the proportion of each stock
        allocate = x_optimal ./ sum(x_optimal); % 1*20
        % Allocate the cash optimal with the proportion
        allocated_cash = cash_optimal * allocate; % 1*20
        % Subtract the cash optimal from original allocated money
        new_money = allocated_money' + allocated_cash; %1*20
        % Find the optimal portion of each stock
        x_optimal = (floor(new_money ./ cur_prices))'; %20*1
        % Calculate the transaction fee
        transaction = cur_prices * abs((x_init - x_optimal)) * 0.005;
        % Calculate the optimal cash account
        cash_optimal = cur_money - cur_prices * x_optimal - transaction;
    end   

end