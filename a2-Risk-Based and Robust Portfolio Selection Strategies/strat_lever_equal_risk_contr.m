function [x_optimal cash_optimal] = strat_lever_equal_risk_contr(x_init, cash_init, mu, Q, cur_prices)
% Leveraged risk parity
% Multiply the number of shares (x_optimal) by 2 only in the first period.
% The money I borrow will be equal to the initial value of portfolio 
% but I will pay interest for each period throughout two years.
    
global cur_year period borrow_money

n = 20;
if cur_year == 2008
    r_rf = 0.045;
else
    r_rf = 0.025;
end
    
cur_money = cur_prices * x_init + cash_init;
period_r = r_rf/6; %interest rate for each period

if period == 1
    borrow_money = cur_money; % borrowed money is equal to initial value of portfolio
    cur_money = cur_money + borrow_money; % take 200% position in equal risk contribution portfolio
else
    cur_money = cur_prices * x_init + cash_init;
    borrow_money = borrow_money;
end

interest = borrow_money * period_r; % borrowed moneny interest for each period

% Equality constraints
A_eq = ones(1,n);
b_eq = 1;

% Inequality constraints
A_ineq = [];
b_ineql = [];
b_inequ = [];
           
% Define initial portfolio
if period == 1
    w0 = 2 * cur_prices' .* x_init / cur_money;
else
    w0 = cur_prices' .* x_init / cur_money;
end

options.lb = zeros(1,n);       % lower bounds on variables
options.lu = ones (1,n);       % upper bounds on variables
options.cl = [b_eq' b_ineql']; % lower bounds on constraints
options.cu = [b_eq' b_inequ']; % upper bounds on constraints

% Set the IPOPT options
options.ipopt.jac_c_constant        = 'yes';
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.mu_strategy           = 'adaptive';
options.ipopt.tol                   = 1e-10;
options.ipopt.print_level           = 0;

% The callback functions
funcs.objective         = @computeObjERC;
funcs.constraints       = @computeConstraints;
funcs.gradient          = @computeGradERC;
funcs.jacobian          = @computeJacobian;
funcs.jacobianstructure = @computeJacobian;

% Run IPOPT
[wsol_2 info] = ipopt(w0',funcs,options);

% Make solution a column vector
if(size(wsol_2,1)==1)
    w_erc = wsol_2';
else
    w_erc = wsol_2;
end

% Risk contribution for the ERC portfolio
% RC_ERC = (w_erc .* ( Q*w_erc )) / sqrt(w_erc'*Q*w_erc);

allocated_money = w_erc * cur_money; % 20*1
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

