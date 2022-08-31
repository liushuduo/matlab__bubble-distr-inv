function [n, varargout] = GlobalInv(f, c_b_inv, bubRadBin, params, Para)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % number of variables
    nBubRadBin = length(bubRadBin);
    problem.nvars = length(bubRadBin) - 1;

    % construct cost function
    freqList = f;                         % all frequency data is known
    % costMetric = @mae;                  % cost metric (mean absolute error)
    % costMetric = @norm;                 % L2 norm 
    costMetric = @(x) norm(x)^2/length(x);         % mean square error

    % randomPerturb = false;               % add random perturb to bubble radius
    randomPerturb = true;
    problem.fitnessfcn = @(x) costFunc(x, c_b_inv, costMetric, bubRadBin, freqList, 1, params, Para, randomPerturb);

    % Linear inequality constraints: sum of all bubble bins less than maximum threshold
    maxBubNum = 1e4;        % maximum bubble number
    problem.Aineq = ones(nBubRadBin-1);  problem.Bineq = ones(nBubRadBin-1, 1) * maxBubNum;

    % Linear equality constraints: none
    problem.Aeq = [];  problem.Beq = [];

    % Lower bound and upper bound;
    problem.lb = zeros(nBubRadBin-1, 1); problem.ub = ones(nBubRadBin-1, 1) * maxBubNum;

    % All variables should be int
    % problem.nonlcon = [];   problem.intcon = 1:length(problem.nvars);

    % GA options
    problem.options = optimoptions('ga', 'Display', 'diagnose');
    % problem.options = optimoptions('ga', 'Display', 'diagnose', 'FunctionTolerance', 1e-7);
    
    tic
    [n, fval, exitflag, output, population, scores] = ga(problem);
    elapsedTime = toc;

    varargout{1} = elapsedTime;
    varargout{2} = fval;
    varargout{3} = exitflag;
    varargout{4} = output;
    varargout{5} = population;
    varargout{6} = scores;

    % save all variables
    save(['global_optimization_ts-', num2str(round(now)), '.mat']);

    % Genetic algorithm fitness function
    function cost = costFunc(bubRadHist, c_b_target, costMetric, bubRadBin, f, regionVol, params, Para, randomPerturb)
        
        % recover bubble radius population from histogram count
        bubRadPop = Radhistc2pop(bubRadHist, bubRadBin, randomPerturb);   
    
        % forward propagation
        c_b = SoundspeedResponse(bubRadPop, f, regionVol, params, Para);
        
        % cost
        cost = costMetric(c_b - c_b_target);
        
        % display cost
    %     disp(['cost: ', num2str(cost)])
    
    end

end

