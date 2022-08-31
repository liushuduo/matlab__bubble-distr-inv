% Definition of physical parameters

% Temperature [Celsius]
Para.T=8.141;

% Gas (Current options: 'CO2' or 'carbon dioxide', 'CH4' or 'methane')
Para.gas='air';

% Current profile
Para.c_slope=0.00;       % Current is zero at z=0, increases linearly at a rate of c_slope m/s per metre, current only in x-direction

% Kinematic viscocity  [ m^2 / s ]
Para.nu=1.4e-6;          % http://www.kayelaby.npl.co.uk/general_physics/2_7/2_7_9.html (average of values at 0 and 20 degrees C)

% Dynamic viscocity [Pa s]
Para.sigma_D=1.5e-3;     % http://www.kayelaby.npl.co.uk/general_physics/2_7/2_7_9.html (average of values at 0 and 20 degrees C)

% Bulk viscocity water [Pa s]
Para.sigma_B=Para.sigma_D*2.2;   % Principles of Sonar Performance Modelling (Mike Ainslie), equation 4.20

% Accr. due to gravity [ m / s^2 ]
Para.g=9.81;

% Seawater density [ kg / m^3 ]
% Para.rhow=1025;
Para.rhow = 1010;

% Atmospheric pressure [Pa]
Para.Patm=101e3;

% Speed of sound in water
% Para.cw=1500;            % [m / s]
Para.cw = 1650;

% Gas constant [m^3 Pa / (K mol) ]
Para.R=8.3144621;

% Surface tension [ N / m ]  % http://web.mit.edu/seawater/Seawater_Property_Tables.pdf (10 degrees, Salinity of 35 ppt)
Para.gamma_st=0.075;        % Equivalent to 75 dyn/cm

switch lower(Para.gas)
    case {'co2','carbon dioxide'}
        % Henry's constant for CO2 in water [ m^3 Pa / mol ]
        Para.H=3100;            % Table 1 c, give H=3.1e4 in atm cm^3 / mol = 3.1e3 in Pa m^3 / mol
        
        % Molar mass of CO2
        Para.Mgas=44/1000;       %  Molar mass of Carbon dioxide [ kg/mol ]
        
        % Speed of sound in CO2
        Para.c_gas=sqrt(4/3*Para.R*(273.15 +Para.T)/Para.Mgas);           % [ m / s ] checked with http://www.phy.mtu.edu/~suits/SpeedofSoundOther.html at 20 degrees C at 1 atm
        
        % Diffusivity (diffusion coefficient) of CO2 in Water [ m^2 / s]
        Para.Df=1.5e-9;          % Table 1 Leifer and Patro, give Df=1.5e-5 in cm^2/s = 1.5e-9 in m^2/s

        % Concentration of CO2 in seawater [ mol / m^3]
%         Para.C=2.3;             % Based on 2.3 mmol/kg taken from http://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Publications/ZeebeWolfEnclp07.pdf
        Para.C=16e-3;            % 10-40 mu mol/kg taken from ''Dissolved Gases and Air-Sea Exchange - ppt'',
                                  % 2e-5 at 0.C, 9.3e-6 at 24 0.C

    case {'ch4','methane'}
                % Henry's constant for CH4 in water [ m^3 Pa / mol ]
        Para.H=79000;            % Table 1 Leifer and Patro, give H=7.9e5 in atm cm^3 / mol = 7.9e4 in Pa m^3 / mol
        
        % Molar mass of CH4
        Para.Mgas=16/1000;       %  Molar mass of Carbon dioxide [ kg/mol ]
        
        % Speed of sound in CO2
        Para.c_gas=sqrt(1.4*Para.R*(273+Para.T)/Para.Mgas);           % [ m / s ] checked with http://www.phy.mtu.edu/~suits/SpeedofSoundOther.html at 20 degrees C at 1 atm

        % Diffusivity (diffusion coefficient) of CH4 in Water [ m^2 / s ]
        Para.Df=1.5e-9;          % Table 1 Leifer and Patro, give Df=1.5e-5 in cm^2/s = 1.5e-9 in m^2/s (same as CO2)
        
        % Concentration of CH4 in seawater [ mol / m^3]
        Para.C=2.5e-6;           % Based on 2.49 n mol/l taken from  http://svr4.terrapub.co.jp/journals/JO/pdf/5101/51010039.pdf (based in Western North Pacific)
    
    case {'o2', 'oxygen'}
        % Henry's constant for O2 in water [ m^3 Pa / mol ]
        Para.H = 76923;         % R. Sander: Compilation of Henry's law constants (version 4.0) for water as solvent, Atmos. Chem. Phys., 15, 4399-4981 (2015), doi:10.5194/acp-15-4399-2015
        
        % Molar mass of O2
        Para.Mgas = 32/1000;    % Molar mass of oxygen [ kg / mol ]

        % Adiabatic constant for oxygen
        Para.gamma = 1.4;

        % Speed of sound in O2
        Para.c_gas = sqrt(Para.gamma*Para.R*(273+Para.T)/Para.Mgas);

        % Diffusivity of oxygen in water [ m^2 / s ]
        Para.Df = 1.9e-9;       % https://www.engineeringtoolbox.com/diffusion-coefficients-d_1404.html

        % Concentration of O2 in seawater [ mol / m^3 ]
        Para.C = 260e-3;        % https://www.soest.hawaii.edu/oceanography/courses/OCN623/Spring2017/12-Non_CO2_gases-2017-web_handouts.pdf

    case {'n2', 'nitrogen'}
        % Henry's constant for N2 in water [ m^3 Pa / mol ]
        Para.H = 156250;    % R. Sander: Compilation of Henry's law constants (version 4.0) for water as solvent, Atmos. Chem. Phys., 15, 4399-4981 (2015), doi:10.5194/acp-15-4399-2015
        
        % Molar mass of N2
        Para.Mgas = 28;     % Molar mass of nitrogen [ kg / mol ]
        
        % Adiabatic constant for nitrogen
        Para.gamma = 1.4;

        % Speed of sound in N2
        Para.c_gas = sqrt(Para.gamma*Para.R*(273+Para.T)/Para.Mgas);
        
        % Diffusivity of nitrogen in water [ m^2 / s ]
        Para.Df = 1.88e-9;  % https://en.wikipedia.org/wiki/Mass_diffusivity
        
        % Concentration of N2 in seawater [ mol / m^3 ]
        Para.C = 470e-3;    % https://www.soest.hawaii.edu/oceanography/courses/OCN623/Spring2017/12-Non_CO2_gases-2017-web_handouts.pdf


    case {'air'}
        
        % Henry's constant for air in water 
        Para.H = 128380;    % 1 / (1/76823 * 21% + 1/15620 * 79%) - weighted averaged by nitrogen and oxygen 
        
        % Molar mass of air
        Para.Mgas = 28.97;

        % Adiabatic constant for nitrogen
        Para.gamma = 1.4;

        % Speed of sound in air
        Para.c_gas = sqrt(Para.gamma*Para.R*(273+Para.T)/Para.Mgas);

        % Diffusivity of air in water [ m^2 / s ]
        Para.Df = 2e-9;
        
        % Concentration of air in seawater [ mol / m^3 ]
        Para.C = 425.9e-3;  % weighted average of oxygen and nitrogen

    otherwise
        error('Unknown gas')
        
end
