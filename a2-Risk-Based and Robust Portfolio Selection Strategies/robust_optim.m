clear all;
clc;
format long;

%% Add PATH to CPLEX and other solvers
addpath('D:\CPLEX\CPLEX1263_x64\cplex\matlab\x64_win64\');

% Random data for 10 stocks
n = 10;
if(exist('portf_data.mat','file')) % read data from file if it exists
    load('portf_data', 'Q', 'mu')
else
    Q = randn(n); Q = Q*Q'/1000; % covariance matrix
    mu  = rand(1,n)'/100;        % expected return
end

%% Define initial portfolio ("equally weighted" or "1/n portfolio")
w0 = ones(n,1) ./ n;

%% Sanity check
Sum_w = sum(w0);
disp(' ')
disp(['For original (equally weighted) portfolio:'])
disp(['  sum of asset weights = ' num2str(Sum_w)]);

%% 1/n portfolio return
ret_init = dot(mu, w0);
disp(['  portfolio return  = ' num2str(ret_init, '%12.10f')])

%% 1/n portfolio variance
var_init = w0' * Q * w0;
disp(['  portfolio st. dev = ' num2str(sqrt(var_init), '%12.10f')])

% Bounds on variables
lb_rMV = zeros(n,1);
ub_rMV = inf*ones(n,1);

% Required portfolio robustness
var_matr = diag(diag(Q));
% Target portfolio return estimation error is return estimation error of 1/n portfolio
rob_init = w0' * var_matr * w0; % return estimation error of initial portfolio
rob_bnd = rob_init; % target return estimation error

% Compute minimum variance portfolio
cplex_minVar = Cplex('MinVar');
cplex_minVar.addCols(zeros(1,n)', [], lb_rMV, ub_rMV);
cplex_minVar.addRows(1, ones(1,n), 1);
cplex_minVar.Model.Q = 2*Q;
cplex_minVar.Param.qpmethod.Cur = 6;
cplex_minVar.solve();
cplex_minVar.Solution
w_minVar = cplex_minVar.Solution.x; % asset weights
ret_minVar = dot(mu, w_minVar);
var_minVar = w_minVar' * Q * w_minVar;
rob_minVar = w_minVar' * var_matr * w_minVar;
disp(['    Return of  minVar portfolio = ' num2str(ret_minVar,'%10.8f')])
disp(['    Return of initial portfolio = ' num2str(ret_init,'%10.8f')])
disp(['  Variance of  minVar portfolio = ' num2str(var_minVar,'%10.8f')])
disp(['  Variance of initial portfolio = ' num2str(var_init,'%10.8f')])
disp(['R.est.err. of  minVar portfolio = ' num2str(sqrt(rob_minVar),'%10.8f')])
disp(['R.est.err. of initial portfolio = ' num2str(sqrt(rob_init),'%10.8f')])
disp(' ')
disp(['Bound on return estimation error of portfolio = ' num2str(rob_bnd,'%10.8f')])
disp(' ')

% Target portfolio return is return of minimum variance portfolio
Portf_Retn = ret_minVar;

%% Formulate and solve robust mean-variance problem
disp(' ')    
disp('======= Solving robust Mean-Variance optimization problem =======')
disp(' ')   
% Objective function
f_rMV  = zeros(n,1);
% Constraints
A_rMV  = sparse([  mu';...
                 ones(1,n)]);
lhs_rMV = [Portf_Retn; 1];
rhs_rMV = [inf; 1];
% Initialize CPLEX environment
cplex_rMV = Cplex('Robust_MV');
% Add objective function and variable bounds
cplex_rMV.addCols(f_rMV, [], lb_rMV, ub_rMV);
% Add constraints
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
cplex_rMV.Solution
    
if(isfield(cplex_rMV.Solution, 'x'))
    w_rMV = cplex_rMV.Solution.x;
    card_rMV = nnz(w_rMV);
    ret_rMV  = dot(mu, w_rMV);
    var_rMV = w_rMV' * Q * w_rMV;
    rob_rMV = w_rMV' * var_matr * w_rMV;
end    
       
disp(' ')
disp(['            Solution status = ' cplex_rMV.Solution.status])
disp(['     Solution status string = ' cplex_rMV.Solution.statusstring])
disp(['              Solution time = ' num2str(cplex_rMV.Solution.time,'%8.3f') ' seconds'])
if(isfield(cplex_rMV.Solution, 'x'))
    disp(['      Portfolio cardinality = ' int2str(card_rMV)])
    disp(['         Solution objective = ' num2str(cplex_rMV.Solution.objval,'%10.8f')])
    disp(['       Sum of asset weights = ' num2str(sum(w_rMV))])
    disp(['       Portfolio rMV return = ' num2str(ret_rMV,'%9.5f')])
    disp(['    Portfolio minVar return = ' num2str(ret_minVar,'%9.5f')])
    disp(['      Portfolio init return = ' num2str(ret_init,'%9.5f')])   
    disp(['      Portfolio rMV st.dev. = ' num2str(sqrt(var_rMV),'%9.5f')])
    disp(['   Portfolio minVar st.dev. = ' num2str(sqrt(var_minVar),'%9.5f')])
    disp(['     Portfolio init st.dev. = ' num2str(sqrt(var_init),'%9.5f')])
    disp(['   Portfolio rMV r.est.err. = ' num2str(sqrt(rob_rMV),'%9.5f')])
    disp([' Portfolio r.est.err. bound = ' num2str(sqrt(rob_bnd),'%9.5f')])
    disp(['Portfolio minVar r.est.err. = ' num2str(sqrt(rob_minVar),'%9.5f')])
    disp(['  Portfolio init r.est.err. = ' num2str(sqrt(rob_init),'%9.5f')])
    disp(['   # of zero portf. weights = ' int2str(sum((w_rMV == 0)))])
    disp(['  # of positive pf. weights = ' int2str(sum((w_rMV > 0)))])
end
disp(' ')

% Round near-zero portfolio weights
w_rMV_nonrnd = w_rMV;
w_rMV(find(w_rMV<=1e-6)) = 0;
w_rMV = w_rMV / sum(w_rMV);
fprintf('\n\nPortfolio weights before and after rounding of near-zero elements:\n')
[w_rMV_nonrnd w_rMV]