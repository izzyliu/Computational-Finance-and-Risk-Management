function  [x_optimal, cash_optimal] = strat_max_Sharpe(x_init, cash_init, mu, Q, cur_prices)
    addpath('/Applications/CPLEX_Studio128/cplex/matlab/x86-64_osx');
    
    % number of stocks
    n = 20;
    
    % Optimization problem data
    lb = zeros(n+1, 1);
    ub = inf*ones(n+1,1);
    
    % Annual risk-free rate for years 2015-2016 is 2.5%
    r_rf = 0.025;
    % Daily risk-free rate
    d_rf = r_rf/252;
    % Add an additional row and an addition column with zeros to Q
    Q = [Q,zeros(20,1);zeros(1,21)];
    % Find the difference between expected return rate and daily risk free rate
    A = [(mu' - d_rf), 0; ones(1,20), -1];
    % Find the bound
    b = [1; 0];
    
    % Compute minimum variance portfolio
    cplex2 = Cplex('max_Sharpe');
    cplex2.addCols(zeros(n+1,1), [], lb, ub);
    cplex2.addRows(b, A, b);
    
    cplex2.Model.Q = 2*Q;
    cplex2.Param.qpmethod.Cur = 6; % concurrent algorithm
    cplex2.Param.barrier.crossover.Cur = 1; % enable crossover
    cplex2.DisplayFunc = []; % disable output to screen
    cplex2.solve();
    
    % Display maximum sharpe ratio portfolio
    sol = cplex2.Solution.x; % 21*1
    y = sol(1:n); % 20*1
    k = sol(21); % Kappa  scalar
    % max sharpe ratio weight for each stock
    w_maxSharpe = y/k; % 20*1
    % Find the total money holding
    cur_money = cur_prices * x_init + cash_init; % scalar
    allocated_money = w_maxSharpe * cur_money; % 20*1
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