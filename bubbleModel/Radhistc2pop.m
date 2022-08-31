function [bubRadPop] = Radhistc2pop(bubRadHist, bubRadBin, randomPerturb)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % recover bubRadHist(bubble radius vector) to bubble radius population
    bubRadInterval = bubRadBin(2:end) - bubRadBin(1:end-1);
    bubRadPop = [];

    for ii = 1 : length(bubRadHist)
        
        n = round(bubRadHist(ii));

        if randomPerturb
            bubRadPop = [bubRadPop, ...
                ones(1, n) * bubRadBin(ii) + rand(1, n) * bubRadInterval(ii)];
        else
            bubRadPop = [bubRadPop, ...
                ones(1, n) * bubRadBin(ii) + linspace(0, bubRadInterval(ii), n)];
        end

    end
    

end