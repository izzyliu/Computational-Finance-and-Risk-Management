function [callMC_Barrier_Knockin_Price_step, putMC_Barrier_Knockin_Price_step] = ...
     MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths)
    
    paths = zeros(numSteps+1, numPaths);
    dT = T/numSteps;
    paths(1,:) = S0;
    
    for iPath = 1:numPaths
        for iStep = 1:numSteps
               paths(iStep+1, iPath) = paths(iStep, iPath) * exp((mu - 0.5*sigma^2)*dT + sigma*sqrt(dT)*normrnd(0,1));                   
        end
    end
    
    knock = sum(paths>= Sb, 1);
    put = zeros(numPaths,1);
    call = zeros(numPaths,1);
    
    for i = 1:numPaths
        if knock(i) > 0
            PutPayoffT = max(K - paths(numSteps+1,i), 0);
            CallPayoffT = max(paths(numSteps+1,i) - K, 0);
            put(i,:) = PutPayoffT * exp(-r*T);
            call(i,:) = CallPayoffT * exp(-r*T);
        else
            call(i,:) = 0;
            put(i,:) = 0;
        end
    end

    callMC_Barrier_Knockin_Price_step = mean(call);
    putMC_Barrier_Knockin_Price_step = mean(put);

end