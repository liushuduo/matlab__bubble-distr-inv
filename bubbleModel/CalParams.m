function [params] = CalParams(freq, params, Para)
%CALPARAMS Summary of this function goes here
%   Detailed explanation goes here
% Calculating general parameters

    [params.ambient.sound_speed,...             % speed of sound [m/s]
        ~,...                                   % dummy
        params.ambient.bulk_modulus,...         % bulk modulus [Pa]
        params.ambient.density,...              % ambient density [kg/m3]
        params.ambient.vapor_pressure,...       % vaport pressure [Pa]
        params.ambient.shear_viscosity,...      % shear viscosity [Pa.s]
        params.ambient.surface_tension] = ...   % surface tension [N/m]
    waterprops( Para.T, Para.S, Para.D, freq );
        
    % Calculating ambient pressure
    params.ambient.pressure = ...               % ambient pressure [Pa]
        Para.Patm + ...
        Para.D*Para.g*Para.rhow;
    
    % Calculating gas density inside bubble (Clay + Medwin 1977 p466)
    params.gas.density = 1.19;%*ones(length(R0),1);
    
    % Calculating thermal diffusivities of gas and ambient liquid
    params.gas.thermal_diffusivity = ...        % gas thermal diffusivity [m2/s]
        params.gas.thermal_conduct./(params.gas.density*params.gas.Cpg);

    params.ambient.thermal_diffusivity = ...    % water thermal diffusivity [m2/s]
        params.ambient.thermal_conduct./(Para.rhow*params.ambient.Cpg);

end

