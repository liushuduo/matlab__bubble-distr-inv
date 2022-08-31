function [n_r, n_i] = ConsLSE(f, c_b, alpha_b, bubRadList, params, Para)
%CONSLSE Summary of this function goes here
%   Detailed explanation goes here

    % recover real part of complex sound speed 
    c_r = c_b;                      

    % recover imaginary part of complex sound speed
    c_i = log(10) / 20 * alpha_b .* c_r.^2 ./ (2*pi * f);

    % the complex sound speed
    c_eff = c_r + 1i * c_i;

    % construct inversion matrix
    K = DeltaKMat(bubRadList, f, params, Para);
    b = (c_eff.^(-2) - Para.cw^(-2)) / Para.rhow;

    K_r = real(K);      K_i = imag(K);
    b_r = real(b);      b_i = imag(b);

    % construct problem
    problem.solver = 'lsqlin';

    % Linear inequality constraints: sum of all bubble bins less than maximum threshold
    maxBubNum = 1e4;        % maximum bubble number
    nBubRadBin = length(bubRadList) + 1;

    % Lower bound and upper bound;
    problem.lb = zeros(nBubRadBin-1, 1); problem.ub = ones(nBubRadBin-1, 1) * maxBubNum;
    
    % problem options
    % problem.options = optimoptions("lsqlin", Algorithm="interior-point", Display="iter", StepTolerance=1e-20);
    problem.options = optimoptions("lsqlin", Algorithm="trust-region-reflective");
    % problem.options = optimoptions('lsqlin', Algorithm='active-set');
    problem.x0 = zeros(nBubRadBin-1, 1);

    % run lsqlin
    % use real/imaginary part to perform least square
    problem.C = real(K);    problem.d = real(b);
    [n_r,resnorm,residual,exitflag,output,lambda] = lsqlin(problem);
    
    problem.C = imag(K);    problem.d = imag(b);
    [n_i,resnorm,residual,exitflag,output,lambda] = lsqlin(problem);

end

