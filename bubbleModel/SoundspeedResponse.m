function [c_b_f, alpha_b_f, c_eff_f] = SoundspeedResponse(bubRadPop, freqList, regionVol, params, Para)
% SOUNDSPEEDRESPONSE Get the sound speed response to bubbly medium

% nulling variables
c_b_f = zeros(size(freqList));            % phase velocity, real part of effective sound speed
alpha_b_f = zeros(size(freqList));        % attenuation, imaginary part of effective sound speed
c_eff_f = zeros(size(freqList));          % effective sound speed

% looping over all frequencies
for iFreq = 1 : length(freqList)
    freq = freqList( iFreq );

    % Calculating general parameters
    params = CalParams(freq, params, Para);

    % nulling var
    vg = 0;
    Delta_K = 0;

    % looping over all bubbles
    for r = bubRadPop
        
        % total gas volume
        vg = vg + 4/3 * pi * r^3;       
    
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
    
        % compute efficient Delta_K
        dtot_res = 2 * dtot_res / ( 2*pi * fr); % M. Hall (1989)
        Delta_K_f_h = 1 / (Para.rhow*pi*(freq^2)) ...
                      * r ...
                      / ((fr/freq)^2 -1 + 1i*dtot_res);
        Delta_K = Delta_K + Delta_K_f_h; % compressibility of bubbles, M. Hall (1989)
    
    end
    
    Delta_K = Delta_K * 1^3 / regionVol;    % convert to delta_K per unit volume
    K_0 = 1 / (Para.rhow * Para.cw^2);      % compressibility of bubbly free water 
    vf = vg / regionVol;                    % void fraction of water column
    K = (1 - vf) * K_0 + Delta_K;           % effective compressibility
    
    % effective density
    params.effective.density = (1-vf)*Para.rhow + vf * params.gas.density;
    C_eff_m2 = K * params.effective.density; 
    c_eff = sqrt( 1 / C_eff_m2 );                   % effective complex sound speed
    c_b = real(c_eff);                              % phase speed
    alpha_b = 20*2*pi*freq/log(10)*imag(1/c_eff);   % Attenuation

    % output data
    c_b_f(iFreq) = c_b;
    alpha_b_f(iFreq) = alpha_b;
    c_eff_f(iFreq) = c_eff;

end

end