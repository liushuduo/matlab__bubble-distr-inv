function [n_r, n_i] = LinearInv(f, c_b, alpha_b, bubRadList, params, Para, tikhonov)
%LINEARINV Summary of this function goes here
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

    n_r = (K_r'*K_r + tikhonov * eye(size(K_r'*K_r))) \ (K_r' * b_r(:));
    n_i = (K_i'*K_i + tikhonov * eye(size(K_i'*K_i))) \ (K_i' * b_i(:));
    
end

