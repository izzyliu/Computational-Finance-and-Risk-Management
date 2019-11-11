function [callMC_European_Price_multi_step, putMC_European_Price_multi_step, paths] = MC_european_price(S0, K, T, r, mu, sigma, numSteps, numPaths)
    % Computes numPaths random paths for a geometric random walk
    % mu is the annual drift, sigma the annual volatility
    % T is the total length of time for the path (in years)
    % dT is the time increment (in years)
    paths = zeros(numSteps+1, numPaths);
    dT = T/numSteps;
    
    % Vector of paths will store realizations of the asset price
    % First asset price is the initial price
    paths(1,:) = S0;
    
    % Generate paths
    for iPath = 1:numPaths
        for iStep = 1:numSteps
               paths(iStep+1, iPath) = paths(iStep, iPath) * exp((mu - 0.5*sigma^2)*dT + sigma*sqrt(dT)*normrnd(0,1));
        end
    end
    
    put = zeros(numPaths,1);
    call = zeros(numPaths,1);
    
    for i = 1:numPaths
        % Calculate the payoff for each path for a Put
        PutPayoffT = max(K - paths(numSteps+1,i), 0);
        % Calculate the payoff for each path for a Call
        CallPayoffT = max(paths(numSteps+1,i) - K, 0);
        % Discount back
        put(i,:) = PutPayoffT * exp(-r*T);
        call(i,:) = CallPayoffT * exp(-r*T);
        
    end
    
    callMC_European_Price_multi_step = mean(call);
    putMC_European_Price_multi_step = mean(put);
    
end