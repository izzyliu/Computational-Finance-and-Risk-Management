clear all;
clc
format long;

Nout  = 100000; % number of out-of-sample scenarios
Nin   = 5000;   % number of in-sample scenarios
Ns    = 5;      % number of idiosyncratic scenarios for each systemic

C = 8;          % number of credit states

% Filename to save out-of-sample scenarios
filename_save_out  = 'scen_out';

% Read and parse instrument data
instr_data = dlmread('instrum_data.csv', ',');
instr_id   = instr_data(:,1);           % ID
driver     = instr_data(:,2);           % credit driver
beta       = instr_data(:,3);           % beta (sensitivity to credit driver)
recov_rate = instr_data(:,4);           % expected recovery rate
value      = instr_data(:,5);           % value
prob       = instr_data(:,6:6+C-1);     % credit-state migration probabilities (default to A)
exposure   = instr_data(:,6+C:6+2*C-1); % credit-state migration exposures (default to A)
retn       = instr_data(:,6+2*C);       % market returns

K = size(instr_data, 1); % number of counterparties

% Read matrix of correlations for credit drivers
rho = dlmread('credit_driver_corr.csv', '\t');
sqrt_rho = (chol(rho))'; % Cholesky decomp of rho (for generating correlated Normal random numbers)

disp('======= Credit Risk Model with Credit-State Migrations =======')
disp('============== Monte Carlo Scenario Generation ===============')
disp(' ')
disp(' ')
disp([' Number of out-of-sample Monte Carlo scenarios = ' int2str(Nout)])
disp([' Number of in-sample Monte Carlo scenarios = ' int2str(Nin)])
disp([' Number of counterparties = ' int2str(K)])
disp(' ')

% Find credit-state for each counterparty
% 8 = AAA, 7 = AA, 6 = A, 5 = BBB, 4 = BB, 3 = B, 2 = CCC, 1 = default
[Ltemp, CS] = max(prob, [], 2); %???????????credit state
clear Ltemp

% Account for default recoveries
exposure(:, 1) = (1-recov_rate) .* exposure(:, 1);% true loss

% Compute credit-state boundaries
CS_Bdry = norminv( cumsum(prob(:,1:C-1), 2) );

% -------- Insert your code here -------- %
% define the number of credit drivers
N_CD = size(rho,1); % 50

if(~exist('scenarios_out.mat','file'))
    
    % -------- Insert your code here -------- %
    % create a matrix y for 100000 scenarios and 50 drivers
    y = zeros(Nout, N_CD); %100000*50
    % create a creditworthiness index matrix w for 100000 scenarios and 100 counterparties
    w = zeros(Nout, K);
    % create a matrix to put the credit state for next year for 100000
    % scenarios and 100 counterparties
    next_CS = zeros(Nout, K);
    % create an idiosyncratic matrix z for 100 counterparties
    z = randn(K,1); % generate a random matrix
    % create a loss matrix for 100000 scenario and 100 counterparties
    Losses_out = zeros(Nout, K);

    for s = 1:Nout
        % -------- Insert your code here -------- %
        random_std_normal = randn(N_CD, 1);
        % create 50 y's for each scenario, since there are 50 credit drivers
        y(s,:) = (sqrt_rho * random_std_normal)';
        for k = 1:K % for each counterparty
            CD = driver(k); % find the corresponding credit driver
            % compute the creditworthiness
            w(s,k) = beta(k) * y(s, CD) + sqrt(1-beta(k)^2)*z(k);
            % put the creditworthiness for scenario s and counterparty k
            % together with the credit state for counterparty k
            % then sorting
            CS_Bdry_w = sort([w(s,k), CS_Bdry(k,:)]);
            CS_index = find(CS_Bdry_w == w(s,k));
            next_CS(s,k) = CS_index;
            % find the corresponding loss for counterparty k and credit
            % state for this counterparty
            Losses_out(s,k) = exposure(k, CS_index);
        end        
    end

    % Calculated out-of-sample losses (100000 x 100)
    % Losses_out

    save('scenarios_out', 'Losses_out')
else
    load('scenarios_out', 'Losses_out')
end

% Normal approximation computed from out-of-sample scenarios
mu_l = mean(Losses_out)'; % compute the mean for each counterparty in 1000000 scenarios 
var_l = cov(Losses_out); % compute the variance for each ounterparty in 100000 scenarios

% Compute portfolio weights
% (1) equal value (dollar amount) is invested in each of 100 bonds
% (2) one unit invested in each of 100 bonds
portf_v = sum(value);     % portfolio value
w0{1} = value / portf_v;  % asset weights (portfolio 1)
w0{2} = ones(K, 1) / K;   % asset weights (portfolio 2)
x0{1} = (portf_v ./ value) .* w0{1};  % asset units (portfolio 1)
x0{2} = (portf_v ./ value) .* w0{2};  % asset units (portfolio 2)

% Quantile levels (99%, 99.9%)
alphas = [0.99 0.999];

% compute the total loss for 100000 scenarios with port1 and port2
total_losses{1} = sort(Losses_out*x0{1});
total_losses{2} = sort(Losses_out*x0{2});
% Compute VaR and CVaR (non-Normal and Normal) for 100000 scenarios
for(portN = 1:2)
    for(q=1:length(alphas))
        alf = alphas(q);
        % -------- Insert your code here -------- %
        VaRout(portN,q) = total_losses{portN}(ceil(Nout*alf));
        VaRinN(portN,q)  = mean(total_losses{portN})+norminv(alf,0,1)*std(total_losses{portN});
        CVaRout(portN,q) = (1/(Nout*(1-alf))) * ( (ceil(Nout*alf)-Nout*alf) * VaRout(portN,q) + sum(total_losses{portN}(ceil(Nout*alf)+1:Nout)) );
        CVaRinN(portN,q) = mean(total_losses{portN}) + (normpdf(norminv(alf,0,1))/(1-alf))*std(total_losses{portN});
        % -------- Insert your code here -------- %        
 end
end


% Perform 100 trials
N_trials = 100;

for(tr=1:N_trials)
    
    % Monte Carlo approximation 1

    % -------- Insert your code here -------- %
    t = 0;
    for s = 1:ceil(Nin/Ns) % systemic scenarios 1000
        % -------- Insert your code here -------- %
        random_std_normal = randn(N_CD, 1);
        % create 50 y's for each scenario, since there are 50 credit drivers
        y_MC1(s,:) = (sqrt_rho * random_std_normal)';
        for si = 1:Ns % idiosyncratic scenarios for each systemic 5
            % -------- Insert your code here -------- %
            % generate 100 z's for each counterparty in each scenario
            z_MC1{s,si} = randn(K,1); 
            for k = 1:K % for each counterparty
                CD = driver(k); % find the corresponding credit driver
                % compute the creditworthiness
                w_MC1{s,si}(k) = beta(k) * y_MC1(s, CD) + sqrt(1-beta(k)^2)*z_MC1{s,si}(k);
                % put the creditworthiness for scenario s and counterparty k
                % together with the credit state for counterparty k
                % then sorting
                CS_Bdry_w_MC1 = sort([w_MC1{s,si}(k), CS_Bdry(k,:)]);
                CS_index_MC1 = find(CS_Bdry_w_MC1 == w_MC1{s,si}(k));
                % find the corresponding loss for counterparty k and credit
                % state for this counterparty
                Losses_out_MC1_temp{s,si}(k) = exposure(k, CS_index_MC1);
            end
            t=t+1;
            Losses_inMC1(t,:) = Losses_out_MC1_temp{s,si};
        end
    end
    
    % Calculated losses for MC1 approximation (5000 x 100)
    % Losses_inMC1
    
    % Monte Carlo approximation 2
    
    % -------- Insert your code here -------- %
    z_MC2 = randn(k,1);
    for s = 1:Nin % systemic scenarios (1 idiosyncratic scenario for each systemic)
        % -------- Insert your code here -------- %
        random_std_normal = randn(N_CD, 1);
        % create 50 y's for each scenario, since there are 50 credit drivers
        y_MC2(s,:) = (sqrt_rho * random_std_normal)';
        for k = 1:K % for each counterparty
                CD = driver(k); % find the corresponding credit driver
                % compute the creditworthiness
                w_MC2(s,k) = beta(k) * y_MC2(s, CD) + sqrt(1-beta(k)^2)*z_MC2(k);
                % put the creditworthiness for scenario s and counterparty k
                % together with the credit state for counterparty k
                % then sorting
                CS_Bdry_w_MC2 = sort([w_MC2(s,k), CS_Bdry(k,:)]);
                CS_index_MC2 = find(CS_Bdry_w_MC2 == w_MC2(s,k));
                % find the corresponding loss for counterparty k and credit
                % state for this counterparty
                Losses_inMC2(s,k) = exposure(k, CS_index_MC2);
        end
    end
        
    % Calculated losses for MC2 approximation (5000 x 100)
    % Losses_inMC2
    
    % Compute VaR and CVaR
    for(portN = 1:2)
        for(q=1:length(alphas))
            alf = alphas(q);
            % -------- Insert your code here -------- %            
            % Compute portfolio loss 
            portf_loss_inMC1{tr,portN} = sort(Losses_inMC1 * x0{portN});
            portf_loss_inMC2{tr,portN} = sort(Losses_inMC2 * x0{portN});
            mu_MC1 = mean(Losses_inMC1)';
            var_MC1 = cov(Losses_inMC1);
            mu_MC2 = mean(Losses_inMC2)';
            var_MC2 = cov(Losses_inMC2);
            % Compute portfolio mean loss mu_p_MC1 and portfolio standard deviation of losses sigma_p_MC1
            mu_p_MC1 = mu_MC1'*x0{portN};
            sigma_p_MC1 = std(portf_loss_inMC1{tr,portN});
            % Compute portfolio mean loss mu_p_MC2 and portfolio standard deviation of losses sigma_p_MC2
            mu_p_MC2 = mu_MC2'*x0{portN};
            sigma_p_MC2 = std(portf_loss_inMC2{tr,portN});
            % Compute VaR and CVaR for the current trial
            VaRinMC1{portN,q}(tr) = portf_loss_inMC1{tr,portN}(ceil(Nin*alf));
            VaRinMC2{portN,q}(tr) = portf_loss_inMC2{tr,portN}(ceil(Nin*alf));
            VaRinN1{portN,q}(tr) = mean(portf_loss_inMC1{tr, portN}) + norminv(alf,0,1)*std(portf_loss_inMC1{tr, portN});
            VaRinN2{portN,q}(tr) = mean(portf_loss_inMC2{tr, portN}) + norminv(alf,0,1)*std(portf_loss_inMC2{tr, portN});
            CVaRinMC1{portN,q}(tr) = (1/(Nin*(1-alf))) * ( (ceil(Nin*alf)-Nin*alf) * VaRinMC1{portN,q}(tr)+ sum(portf_loss_inMC1{tr, portN}(ceil(Nin*alf)+1:Nin)));
            CVaRinMC2{portN,q}(tr) = (1/(Nin*(1-alf))) * ( (ceil(Nin*alf)-Nin*alf) * VaRinMC2{portN,q}(tr)+ sum(portf_loss_inMC2{tr, portN}(ceil(Nin*alf)+1:Nin)));
            CVaRinN1{portN,q}(tr) = mean(portf_loss_inMC1{tr, portN}) + (normpdf(norminv(alf,0,1))/(1-alf))*std(portf_loss_inMC1{tr, portN});
            CVaRinN2{portN,q}(tr) = mean(portf_loss_inMC2{tr, portN}) + (normpdf(norminv(alf,0,1))/(1-alf))*std(portf_loss_inMC2{tr, portN});
            % -------- Insert your code here -------- %
        end
    end
end

% Display portfolio VaR and CVaR
for(portN = 1:2)
fprintf('\nPortfolio %d:\n\n', portN)    
 for(q=1:length(alphas))
    alf = alphas(q);
    fprintf('Out-of-sample: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, VaRout(portN,q), 100*alf, CVaRout(portN,q))
    fprintf('In-sample MC1: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinMC1{portN,q}), 100*alf, mean(CVaRinMC1{portN,q}))
    fprintf('In-sample MC2: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinMC2{portN,q}), 100*alf, mean(CVaRinMC2{portN,q}))
    fprintf(' In-sample No: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, VaRinN(portN,q), 100*alf, CVaRinN(portN,q))
    fprintf(' In-sample N1: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinN1{portN,q}), 100*alf, mean(CVaRinN1{portN,q}))
    fprintf(' In-sample N2: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n\n', 100*alf, mean(VaRinN2{portN,q}), 100*alf, mean(CVaRinN2{portN,q}))
 end
end

% Plot portfolio 1 with out-of-sample
figure(1);
% -------- Insert your code here -------- %
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(total_losses{1}, 100);
bar(binLocations, frequencyCounts);
hold on;
line([VaRout(1,1) VaRout(1,1)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
line([VaRout(1,2) VaRout(1,2)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
line([CVaRout(1,1) CVaRout(1,1)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
line([CVaRout(1,2) CVaRout(1,2)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;

normf = ( 1/(std(total_losses{1})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(total_losses{1}))/std(total_losses{1})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1); 
hold on;

text(0.98*VaRout(1,1), max(frequencyCounts)/1.9, '99% VaR')
text(0.98*VaRout(1,2), max(frequencyCounts)/1.9, '99.9% VaR')
text(0.98*CVaRout(1,1), max(frequencyCounts)/1.9, '99% CVaR')
text(0.98*CVaRout(1,2), max(frequencyCounts)/1.9, '99.9% CVaR')
hold off;
xlabel('Credit-State Migration Losses')
ylabel('Frequency')
title('Out-of-sample Distribution for Portfolio 1')

% Plot portfolio 2 with out-of-sample
figure(2);
% -------- Insert your code here -------- %
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(total_losses{2}, 100);
bar(binLocations, frequencyCounts);
hold on;
line([VaRout(2,1) VaRout(2,1)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
line([VaRout(2,2) VaRout(2,2)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
line([CVaRout(2,1) CVaRout(2,1)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
line([CVaRout(2,2) CVaRout(2,2)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;

normf = ( 1/(std(total_losses{2})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(total_losses{2}))/std(total_losses{2})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1); 
hold on;

text(0.98*VaRout(1,1), max(frequencyCounts)/1.9, '99% VaR')
text(0.98*VaRout(1,2), max(frequencyCounts)/1.9, '99.9% VaR')
text(0.98*CVaRout(1,1), max(frequencyCounts)/1.9, '99% CVaR')
text(0.98*CVaRout(1,2), max(frequencyCounts)/1.9, '99.9% CVaR')
hold off;
xlabel('Credit-State Migration Losses')
ylabel('Frequency')
title('Out-of-sample Distribution for Portfolio 2')

% Plot portfolio 1 with in sample MC1
VaRinMC1_99 = mean(VaRinMC1{1,1});
VaRinMC1_999 = mean(VaRinMC1{1,2});
CVaRinMC1_99 = mean(CVaRinMC1{1,1});
CVaRinMC1_999 = mean(CVaRinMC1{1,2});

% compute the average loss for each scenario
p1_loss_inMC1 = zeros(5000,1);
for i = 1:100
    p1_loss_inMC1 = p1_loss_inMC1 + portf_loss_inMC1{i,1};
end
mean_p1_loss_inMC1 = p1_loss_inMC1.*(1/100);

figure(3);
set(gcf,'color', 'white');
[frequencyCounts, binLocations] = hist(mean_p1_loss_inMC1, 100);
bar(binLocations, frequencyCounts);
hold on;
line([VaRinMC1_99 VaRinMC1_99], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
line([VaRinMC1_999 VaRinMC1_999], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
line([CVaRinMC1_99 CVaRinMC1_99], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
line([CVaRinMC1_999 CVaRinMC1_999], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;

normf = ( 1/(std(mean_p1_loss_inMC1)*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(mean_p1_loss_inMC1))/std(mean_p1_loss_inMC1)).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1);
hold on;

text(0.98*VaRinMC1_99, max(frequencyCounts)/1.9, '99% VaR')
text(0.98*VaRinMC1_999, max(frequencyCounts)/1.9, '99.9% VaR')
text(0.98*CVaRinMC1_99, max(frequencyCounts)/1.9, '99% CVaR')
text(0.98*CVaRinMC1_999, max(frequencyCounts)/1.9, '99.9% CVaR')
hold off;
xlabel('Credit-State Migration Losses')
ylabel('Frequency')
title('Monte Carlo approx. 1 Distribution for Portfolio 1')

% Plot portfolio 2 with in sample MC1
VaRinMC1_99_2 = mean(VaRinMC1{2,1});
VaRinMC1_999_2 = mean(VaRinMC1{2,2});
CVaRinMC1_99_2 = mean(CVaRinMC1{2,1});
CVaRinMC1_999_2 = mean(CVaRinMC1{2,2});

% compute the average loss for each scenario
p2_loss_inMC1 = zeros(5000,1);
for i = 1:100
    p2_loss_inMC1 = p2_loss_inMC1 + portf_loss_inMC1{i,2};
end
mean_p2_loss_inMC1 = p2_loss_inMC1.*(1/100);

figure(4);
set(gcf,'color', 'white');
[frequencyCounts, binLocations] = hist(mean_p2_loss_inMC1, 100);
bar(binLocations, frequencyCounts);
hold on;
line([VaRinMC1_99_2 VaRinMC1_99_2], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
line([VaRinMC1_999_2 VaRinMC1_999_2], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
line([CVaRinMC1_99_2 CVaRinMC1_99_2], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
line([CVaRinMC1_999_2 CVaRinMC1_999_2], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;

normf = ( 1/(std(mean_p2_loss_inMC1)*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(mean_p2_loss_inMC1))/std(mean_p2_loss_inMC1)).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1);
hold on;

text(0.98*VaRinMC1_99_2, max(frequencyCounts)/1.9, '99% VaR')
text(0.98*VaRinMC1_999_2, max(frequencyCounts)/1.9, '99.9% VaR')
text(0.98*CVaRinMC1_99_2, max(frequencyCounts)/1.9, '99% CVaR')
text(0.98*CVaRinMC1_999_2, max(frequencyCounts)/1.9, '99.9% CVaR')
hold off;
xlabel('Credit-State Migration Losses')
ylabel('Frequency')
title('Monte Carlo approx. 1 Distribution for Portfolio 2')

% Plot portfolio 1 with in sample MC2
VaRinMC2_99 = mean(VaRinMC2{1,1});
VaRinMC2_999 = mean(VaRinMC2{1,2});
CVaRinMC2_99 = mean(CVaRinMC2{1,1});
CVaRinMC2_999 = mean(CVaRinMC2{1,2});

% compute the average loss for each scenario
p1_loss_inMC2 = zeros(5000,1);
for i = 1:100
    p1_loss_inMC2 = p1_loss_inMC2 + portf_loss_inMC2{i,1};
end
mean_p1_loss_inMC2 = p1_loss_inMC2.*(1/100);

figure(5);
set(gcf,'color', 'white');
[frequencyCounts, binLocations] = hist(mean_p1_loss_inMC2, 100);
bar(binLocations, frequencyCounts);
hold on;
line([VaRinMC2_99 VaRinMC2_99], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
line([VaRinMC2_999 VaRinMC2_999], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
line([CVaRinMC2_99 CVaRinMC2_99], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
line([CVaRinMC2_999 CVaRinMC2_999], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;

normf = ( 1/(std(mean_p1_loss_inMC2)*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(mean_p1_loss_inMC2))/std(mean_p1_loss_inMC2)).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1);
hold on;

text(0.98*VaRinMC2_99, max(frequencyCounts)/1.9, '99% VaR')
text(0.98*VaRinMC2_999, max(frequencyCounts)/1.9, '99.9% VaR')
text(0.98*CVaRinMC2_99, max(frequencyCounts)/1.9, '99% CVaR')
text(0.98*CVaRinMC2_999, max(frequencyCounts)/1.9, '99.9% CVaR')
hold off;
xlabel('Credit-State Migration Losses')
ylabel('Frequency')
title('Monte Carlo approx. 2 Distribution for Portfolio 1')

% Plot portfolio 2 with in sample MC2
VaRinMC2_99_2 = mean(VaRinMC2{2,1});
VaRinMC2_999_2 = mean(VaRinMC2{2,2});
CVaRinMC2_99_2 = mean(CVaRinMC2{2,1});
CVaRinMC2_999_2 = mean(CVaRinMC2{2,2});

% compute the average loss for each scenario
p2_loss_inMC2 = zeros(5000,1);
for i = 1:100
    p2_loss_inMC2 = p2_loss_inMC2 + portf_loss_inMC2{i,2};
end
mean_p2_loss_inMC2 = p2_loss_inMC2.*(1/100);

figure(6);
set(gcf,'color', 'white');
[frequencyCounts, binLocations] = hist(mean_p2_loss_inMC2, 100);
bar(binLocations, frequencyCounts);
hold on;
line([VaRinMC2_99_2 VaRinMC2_99_2], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
line([VaRinMC2_999_2 VaRinMC2_999_2], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
line([CVaRinMC2_99_2 CVaRinMC2_99_2], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
line([CVaRinMC2_999_2 CVaRinMC2_999_2], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;

normf = ( 1/(std(mean_p2_loss_inMC2)*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(mean_p2_loss_inMC2))/std(mean_p2_loss_inMC2)).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1);
hold on;

text(0.98*VaRinMC2_99_2, max(frequencyCounts)/1.9, '99% VaR')
text(0.98*VaRinMC2_999_2, max(frequencyCounts)/1.9, '99.9% VaR')
text(0.98*CVaRinMC2_99_2, max(frequencyCounts)/1.9, '99% CVaR')
text(0.98*CVaRinMC2_999_2, max(frequencyCounts)/1.9, '99.9% CVaR')
hold off;
xlabel('Credit-State Migration Losses')
ylabel('Frequency')
title('Monte Carlo approx. 2 Distribution for Portfolio 2')

