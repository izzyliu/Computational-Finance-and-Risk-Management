function  [x_optimal, cash_optimal] = strat_equally_weighted(x_init, cash_init, mu, Q, cur_prices)

    % find the initial total money holding
    cur_money = cur_prices * x_init + cash_init;
    % allocate the initial money to 20 stocks equally
    equal_w = cur_money * (1/20);
    equal_money = ones(20,1) * equal_w;

    x_optimal = round(equal_money./cur_prices');
    transaction = cur_prices * abs((x_init - x_optimal)) * 0.005;
    cash_optimal = cur_money - cur_prices * x_optimal - transaction;

    % check if the cash account is non_negative.
    if cash_optimal < 0
        % reserve the amount of negative cash account before the investment,
        % subtract from the initial total money.
        new_money = cur_money + cash_optimal;
        new_equal_money = ones(20,1) * new_money * (1/20);
    
        x_optimal = floor(new_equal_money./cur_prices');
        transaction = cur_prices * abs((x_init - x_optimal)) * 0.005;
        cash_optimal = cur_money - cur_prices * x_optimal - transaction;
    end

end