function [leftmoney] = gameround(leftmoney,weight,winrate)
    a = rand();
    if a < winrate
        leftmoney = leftmoney + weight;
    else
        leftmoney = leftmoney - weight;
    end
end

