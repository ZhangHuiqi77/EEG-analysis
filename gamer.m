function [money] = gamer(initialM,initialW,winRatio,maxTurn)
    % Here we want to simulate a game. Give the intial money, initial weight,
    % win ratio and the maximum turn as input. Then, if the gamer win in this
    % round, add the weight to the leftMoney and keep the weight as initial
    % weight. If the gamer lose the round, give 2 times of initial weight as
    % current weight, if the gammer has enough money left. After every round,
    % update the current money in the variable 'money'.
    curMoney = initialM;
    curWeight = initialW;
    money = zeros(maxTurn+1,1);
    money(1) = initialM;
    for m = 1:maxTurn
        leftMoney = gameround(curMoney,curWeight,winRatio);% call functions gameround
        if leftMoney ==0
            money((m+1):maxTurn) = [];
            return;
        end
        if leftMoney > curMoney
            curWeight = initialW;
        else
            if leftMoney > 2*curWeight;
                curWeight = 2*curWeight;
            else
                curWeight = leftMoney;
            end
        end
        curMoney = leftMoney;
        money(m+1) = curMoney;
    end
end

