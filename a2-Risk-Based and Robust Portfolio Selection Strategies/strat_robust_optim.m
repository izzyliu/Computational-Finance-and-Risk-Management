function [x_optimal cash_optimal] = strat_robust_optim(x_init, cash_init, mu, Q, cur_prices)

% Random data for 20 stocks
n = 20;

% total money
cur_money = cur_prices * x_init + cash_init;

% define the initial weight
w0 = cur_prices' .* x_init / cur_money;

% Bounds on variables
lb_rMV = zeros(n,1); 
ub_rMV = inf*ones(n,1);

% Target portfolio return estimation error
var_matr = diag(diag(Q));
rob_init = w0' * var_matr * w0; % return estimate error of 1/n portf 
rob_bnd = rob_init; % target return estimation error

% Compute minimum variance portfolio (MVP)
% Compute minimum variance portfolio
cplex_minVar = Cplex('MinVar');
cplex_minVar.addCols(zeros(1,n)', [], lb_rMV, ub_rMV);
cplex_minVar.addRows(1, ones(1,n), 1);
cplex_minVar.Model.Q = 2*Q;
cplex_minVar.Param.qpmethod.Cur = 6;
cplex_minVar.DisplayFunc = []; % disable output to screen
cplex_minVar.solve();
cplex_minVar.Solution;
w_minVar = cplex_minVar.Solution.x; % asset weights
ret_minVar = dot(mu, w_minVar);
var_minVar = w_minVar' * Q * w_minVar;
rob_minVar = w_minVar' * var_matr * w_minVar;

% Target portfolio return = return of MVP 
Portf_Retn = ret_minVar;

% Formulate and solve robust mean-variance problem
% objective function
f_rMV = zeros(n,1); 
% Constraints
A_rMV = sparse([mu'; ones(1,n)]);
lhs_rMV = [Portf_Retn; 1]; 
rhs_rMV = [inf; 1];
% Create CPLEX model
cplex_rMV = Cplex('Robust_MV');
cplex_rMV.addCols(f_rMV, [], lb_rMV, ub_rMV);
cplex_rMV.addRows(lhs_rMV, A_rMV, rhs_rMV);
% Add quadratic objective
cplex_rMV.Model.Q = 2*Q;
% Add quadratic constraint on return estimation error (robustness constraint)
Qq_rMV = var_matr;
cplex_rMV.addQCs(zeros(size(f_rMV)), Qq_rMV, 'L', rob_bnd, {'qc_robust'});
% Set CPLEX parameters
cplex_rMV.Param.threads.Cur = 4;
cplex_rMV.Param.timelimit.Cur = 60;
cplex_rMV.Param.barrier.qcpconvergetol.Cur = 1e-12; % solution tolerance
cplex_rMV.DisplayFunc = []; % disable output to screen 
cplex_rMV.solve();   
cplex_rMV.Solution;
    
if(isfield(cplex_rMV.Solution, 'x'))
    w_rMV = cplex_rMV.Solution.x;
    card_rMV = nnz(w_rMV);
    ret_rMV  = dot(mu, w_rMV);
    var_rMV = w_rMV' * Q * w_rMV;
    rob_rMV = w_rMV' * var_matr * w_rMV;
end

% Round near-zero portfolio weights
w_rMV_nonrnd = w_rMV;
w_rMV(find(w_rMV<=1e-6)) = 0;
w_rMV = w_rMV / sum(w_rMV);

allocated_money = w_rMV * cur_money; % 20*1
x_optimal = round(allocated_money' ./ cur_prices)'; % 20*1
transaction = cur_prices * abs((x_init - x_optimal)) * 0.005; %scalar
cash_optimal = cur_money - cur_prices * x_optimal - transaction; %scalar
    
if cash_optimal < 0
    % Find the proportion of each stock
    allocate = (x_optimal ./ sum(x_optimal))'; % 1*20
    % Allocate the cash optimal with the proportion
    allocated_cash = cash_optimal * allocate; % 1*20
    % Subtract the cash optimal from original allocated money
    new_money = allocated_money + allocated_cash'; %20*1
    % Find the optimal portion of each stock
    x_optimal = floor(new_money ./ cur_prices'); %20*1
    % Calculate the transaction fee
    transaction = cur_prices * abs((x_init - x_optimal)) * 0.005;
    % Calculate the optimal cash account
    cash_optimal = cur_money - cur_prices * x_optimal - transaction;
end

end

