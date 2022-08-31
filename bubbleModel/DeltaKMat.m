function [A] = DeltaKMat(bubRadList, freqList, params, Para)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Nulling matrix
A = zeros(length(freqList), length(bubRadList));

for iFreq = 1 : length(freqList)
    freq = freqList( iFreq );

    % Calculating general parameters
    params = CalParams(freq, params, Para);
    
    for iRad = 1 : length(bubRadList)
        r = bubRadList( iRad );

        % bubble resonance frequency (Minneart)
        fr = 1/(2*pi*r)* sqrt(3 * params.gas.ratio_spec_heat*params.ambient.pressure/Para.rhow);
        
        % compute damping coefficient
        [dtot_res, ...           % total damping coefficient
        ~, ...	% thermal damping coefficient
        ~, ...	% viscous damping coefficient
        ~, ...	% radiation damping coefficient
        ] = ...          % polytropic exponent
        bubble_damping_properetti(fr,...
        r,...
        Para.cw,...
        Para.c_gas,...
        params.ambient.thermal_diffusivity,...
        params.gas.thermal_diffusivity,...
        params.ambient.thermal_conduct,...
        params.gas.thermal_conduct,...
        params.gas.ratio_spec_heat,...
        Para.rhow,...
        params.gas.density,...
        params.ambient.shear_viscosity);
        
        dtot_res = 2 * dtot_res / ( 2*pi * fr); % M. Hall (1989)

        % compute efficient Delta_K for current bin
        Delta_K_f_h = 1 / (Para.rhow*pi*(freq^2)) ...
                      * r ...
                      / ((fr/freq)^2 -1 + 1i*dtot_res);

        A(iFreq, iRad) = Delta_K_f_h;
    end
end

end