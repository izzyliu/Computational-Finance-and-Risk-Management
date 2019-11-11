clc;
clear all;
format long

% Pricing a European option using Black-Scholes formula and Monte Carlo simulations
% Pricing a Barrier option using Monte Carlo simulations

S0 = 100;     % spot price of the underlying stock today
K = 105;      % strike at expiry
mu = 0.05;    % expected return
sigma = 0.18;  % volatility
r = 0.05;     % risk-free rate
T = 1.0;      % years to expiry
Sb = 110;     % barrier
numPaths = 5000; % number of paths
numSteps = 252;  % number of steps
%numSteps = 12;
%numSteps = 24;

% Define variable numSteps to be the number of steps for multi-step MC
% numPaths - number of sample paths used in simulations

% Implement your Black-Scholes pricing formula
[call_BS_European_Price, putBS_European_Price] = BS_european_price(S0, K, T, r, sigma);

% Implement your one-step Monte Carlo pricing procedure for European option
% numSteps = 1;
[callMC_European_Price_1_step, putMC_European_Price_1_step] = MC_european_price(S0, K, T, r, mu, sigma, 1, numPaths);

% Implement your multi-step Monte Carlo pricing procedure for European option
[callMC_European_Price_multi_step, putMC_European_Price_multi_step, paths] = MC_european_price(S0, K, T, r, mu, sigma, numSteps, numPaths);

% Implement your one-step Monte Carlo pricing procedure for Barrier option
% numSteps = 1;
[callMC_Barrier_Knockin_Price_1_step, putMC_Barrier_Knockin_Price_1_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, 1, numPaths);

% Implement your multi-step Monte Carlo pricing procedure for Barrier option
[callMC_Barrier_Knockin_Price_multi_step, putMC_Barrier_Knockin_Price_multi_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths);

disp(['Black-Scholes price of an European call option is ',num2str(call_BS_European_Price)])
disp(['Black-Scholes price of an European put option is ',num2str(putBS_European_Price)])
disp(['One-step MC price of an European call option is ',num2str(callMC_European_Price_1_step)])
disp(['One-step MC price of an European put option is ',num2str(putMC_European_Price_1_step)])
disp(['Multi-step MC price of an European call option is ',num2str(callMC_European_Price_multi_step)])
disp(['Multi-step MC price of an European put option is ',num2str(putMC_European_Price_multi_step)])
disp(['One-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_1_step)])
disp(['One-step MC price of an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_1_step)])
disp(['Multi-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_multi_step)])
disp(['Multi-step MC price of an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_multi_step)])


% Plot results
figure(1);
%%%%%%%%%%% Insert your code here %%%%%%%%%%%%
%Plot one chart illustrates your Monte Carlo pricing procedure in the best way.
for i=1:numPaths 
    plot(0:numSteps,paths(:,i));
    hold on;
end
hold off;
%axis([0 numSteps+1 30 inf])
title('Monte Carlo pricing procedure for European option');

numPathsList = [100,500,1000,5000,10000,20000,30000,40000,50000];
numStepsList = [2, 12, 24, 252];
step_path = zeros(length(numStepsList), length(numPathsList));

for iPath = 1:length(numPathsList)
    for iStep = 1:length(numStepsList)
        [call, put] = MC_european_price(S0, K, T, r, mu, sigma, numStepsList(iStep), numPathsList(iPath));
        step_path(iStep, iPath) = abs(call-call_BS_European_Price)+abs(put-putBS_European_Price);
    end
end
minimum = min(min(step_path));
[x,y]=find(step_path==minimum);
step = numStepsList(x);
path = numPathsList(y);
disp(['The number of steps is ',num2str(step)])
disp(['The number of paths is ',num2str(path)])
        