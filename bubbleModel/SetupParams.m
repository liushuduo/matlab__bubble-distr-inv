function [params] = SetupParams(gasName, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

p = inputParser();

% general parameters
addParameter(p, 'atmospheric_pressure', 101325)      % atmospheric pressure [Pa]
addParameter(p, 'gravity', 9.80665)              % acceleration due to gravity [m/s2]
addParameter(p, 'free_gas_dens_sealevel', 1.29)      % density of free gas at sea level for air [kg/m^3] @ Clay + Medwin 1977 p466
addParameter(p, 'gas_constant', 8.31441)             % universal gas constant [J/(K*mol)]

% ambient parameters
addParameter(p, 'ambient_salinity', 0)               % ambient salinity [ppt]
addParameter(p, 'ambient_Cpv', 4.193e3);             % water heat capacity at constant volume [J/(kg*K)]
addParameter(p, 'ambient_Cpg', 4.193e3);             % water heat capacity at constant pressure [J/(kg*K)]   http://www.engineeringtoolbox.com/water-thermal-properties-d_162.html
addParameter(p, 'ambient_thermal_conduct', 56.3e-2)  % water thermal conductivity


% gas parameters
switch lower(gasName)

    case {'co2', 'carbon dioxide'}
        gas.Cpg = 0.834e3;                       % gas heat capacity at constant pressure [J/(kg*K)]
        gas.thermal_conduct = 10.5e-2;           % gas thermal conductivity [W/(m*K)]
        gas.mol_weight = 44.01;                  % gas molecular weight [g/mol]
        gas.ratio_spec_heat = 1.304;             % ratio of specific heat of a gas at constant pressure to that at constant volume
        gas.sound_speed = 259;                   % sound speed in CO2 [m/s]
    
    case 'air'
        gas.Cpg = 1.005;                        % air heat capacity at constant pressure [J/(kg*K)]
        gas.thermal_conduct = 25.24e-3;         % air thermal conductivity [W/(m*K)]
        gas.mol_weight = 28.97;                 % air molecular weight [g/mol]
        gas.ratio_spec_heat = 1.4;
        gas.sound_speed = 334;    
    
    otherwise
        error('unknown gas.')
end

% parse parameters
parse(p, varargin{:})
params.atmospheric_pressure = p.Results.atmospheric_pressure;
params.gravity_acc = p.Results.gravity;
params.free_gas_dens_sealevel = p.Results.free_gas_dens_sealevel;
params.gas_constant = p.Results.gas_constant;

params.ambient.salinity = p.Results.ambient_salinity;
params.ambient.Cpv = p.Results.ambient_Cpv;
params.ambient.Cpg = p.Results.ambient_Cpg;
params.ambient.thermal_conduct = p.Results.ambient_thermal_conduct;

params.gas = gas;
end