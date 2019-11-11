function  [x_optimal, cash_optimal] = strat_buy_and_hold_equally(x_init, cash_init, mu, Q, cur_prices)
    % find the initial total money holding
    cur_money = cur_prices * x_init + cash_init;
    x_optimal = [1061;3230;1493;94;2715;989;754;1056;453;306;1465;1797;2772;1364;18587;2413;2465;160;1282;1226]; % 20*1
    transaction = cur_prices * abs((x_init - x_optimal)) * 0.005;
    cash_optimal = cur_money - cur_prices * x_optimal - transaction;

    % check if the cash account is non_negative.
    if cash_optimal < 0
       % Find the proportion of each stock
        allocate = x_optimal' ./ sum(x_optimal); % 1*20
        % Allocate the cash optimal with the proportion
        allocated_cash = cash_optimal * allocate; % 1*20
        % Subtract the cash optimal from original allocated money
        position_change = floor(allocated_cash ./ cur_prices); % 1*20
        % Find the optimal portion of each stock
        x_optimal = x_optimal - position_change'; %20*1
        % Calculate the transaction fee
        transaction = cur_prices * abs((x_init - x_optimal)) * 0.005;
        % Calculate the optimal cash account
        cash_optimal = cur_money - cur_prices * x_optimal - transaction;
    end

end