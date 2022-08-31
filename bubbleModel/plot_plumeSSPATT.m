
% frequency 
% clear all
% close all
% 
% load bubbleCH4plumerate5dt1stime500swide.mat;
% define_phys_parameters;
% load S:\MATLAB\matlab__soundspeed-analysis\bubble-model\zhoushan_flux-rate-0.2_bubpdf-experiment.mat
%%

Para.S = 35;                                    % salinity [ppt]
Para.D = depth;

% general parameters
params.atmospheric_pressure = 101325;           % atmospheric pressure [Pa]
params.gravity_acc = 9.80665;                   % acceleration due to gravity [m/s2]
params.free_gas_dens_sealevel = 1.29;           % density of free gas at sea level for air [kg/m^3] @ Clay + Medwin 1977 p466
params.gas_constant = 8.31441;                  % universal gas constant [J/(K*mol)]

% ambient parameters
params.ambient.salinity = 35;                   % salinity [ppt]
params.ambient.Cpv = 4.193e3;                   % water heat capacity at constant volume [J/(kg*K)]
params.ambient.Cpg = 4.193e3;                   % water heat capacity at constant pressure [J/(kg*K)]   http://www.engineeringtoolbox.com/water-thermal-properties-d_162.html
params.ambient.thermal_conduct = 56.3e-2;       % water thermal conductivity

% gas parameters
params.gas.Cpg = 0.834e3;                       % gas heat capacity at constant pressure [J/(kg*K)]
params.gas.thermal_conduct = 10.5e-2;          % gas thermal conductivity [W/(m*K)]
params.gas.mol_weight = 44.01;                  % gas molecular weight
params.gas.ratio_spec_heat = 1.304;             % ratio of specific heat of a gas at constant pressure to that at constant volume
params.gas.sound_speed = 259;                   % sound speed in CO2 [m/s]

f = 100:100:25000;
c_eff = zeros(length(f),1);
c_b = zeros(length(f),1);
alpha_b = zeros(length(f),1);

for countf = 1:length(f) % acoustic frequency in Hz
    disp(countf);
    ff = f(countf);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Calculating general parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % calculating general parameters
    [params.ambient.sound_speed,...             % speed of sound [m/s]
        ~,...                                   % dummy
        params.ambient.bulk_modulus,...         % bulk modulus [Pa]
        params.ambient.density,...              % ambient density [kg/m3]
        params.ambient.vapor_pressure,...       % vaport pressure [Pa]
        params.ambient.shear_viscosity,...      % shear viscosity [Pa.s]
        params.ambient.surface_tension] = ...   % surface tension [N/m]
        waterprops( Para.T,...
        Para.S,...
        Para.D,...
        ff);
    %     sf
    % calculating ambient pressure
    params.ambient.pressure = ...               % ambient pressure [Pa]
        Para.Patm + ...
        Para.D*Para.g*Para.rhow;
    
    % calculating internal bubble pressure
    % params.bubble.internal_pressure = ...       % bubble internal pressure [Pa]
    %     params.ambient.pressure + 2*params.ambient.surface_tension./R0;
    
    % calculating gas density inside bubble (Clay + Medwin 1977 p466)
    params.gas.density = 1.19;%*ones(length(R0),1);
    
    % calculating thermal diffusivities of gas and ambient liquid
    params.gas.thermal_diffusivity = ...        % gas thermal diffusivity [m2/s]
        params.gas.thermal_conduct./(params.gas.density*params.gas.Cpg);
    
    params.ambient.thermal_diffusivity = ...    % water thermal diffusivity [m2/s]
        params.ambient.thermal_conduct./(Para.rhow*params.ambient.Cpg);
    
    
    Nb=length(Bubbles);
    
    scale=12;
    
    Nbubs=length(Bubbles);
    
    theta=linspace(0,2*pi,100);
    xc=sin(theta);
    yc=cos(theta);
    
    % Tvalues=[1,2,5,10,15,20];
    %     Tvalues=20;
    
    for k=1:Nbubs
        if (isempty(Bubbles(k).tstop)), Bubbles(k).tstop=inf; end
    end
    
    close
    % for count=1:length(Tvalues)
    %     T=Tvalues(count);
    
    T = 50;
    hh = 1;
    Vg = 0;
%     Vf = zeros(length(h),1);
%     kr = [];
    
    % for h = 0.1:0.1:3
    indR0 = 0;
%     R0 = [];
%     Fr = [];
    Delta_K = 0;
    % figure(100)
    
    for k=1:Nbubs
        
        if ( Bubbles(k).tstart<T && T<Bubbles(k).tstop )
            m=round((T-Bubbles(k).tstart)/dt)+1;
            m=max([m 1]);
        else
            continue
        end
        
        r=Bubbles(k).r(m);
        x=Bubbles(k).x(m);
        z=Bubbles(k).z(m);
        %        subplot(1,1,count),plot(xc*r*scale+x*0.25,yc*r*scale+z),hold on,axis equal,axis([-1 1 0 3])
        if ( z <= hh && z > hh - 1 )
            z;
            indR0 = indR0 +1;
%             R0(indR0) = r;
%             kr = [kr k];
            Vg = Vg + 4/3 * pi * r^3; % Gas volume
            r;
            %         size(r)
            %         pause
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Calculating bubble resonances
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % calculating resonant frequencies using leighton model
            fr = 1/(2*pi*r)* sqrt(3*params.gas.ratio_spec_heat*params.ambient.pressure/...
                Para.rhow);         % Minnaert resonant frequencies
            %             fr;
            % calculating resonant frequencies using leighton model
%             fr = ...                         % resonant frequencies
%                 bubble_resonance_leighton(  r,...
%                 params.ambient.shear_viscosity,...
%                 params.gas.ratio_spec_heat,...
%                 params.ambient.vapor_pressure,...
%                 params.ambient.surface_tension,...
%                 params.ambient.pressure,...
%                 params.ambient.density,...
%                 1);
            %             fr;
%             Fr(indR0) = fr;
            %         size(fr)
            %         pause
            %             f0(indR0) = fr;
            %                 params.ambient.sound_speed,...     % two cases, (given, calculated)
            %                 params.gas.sound_speed,...         % two cases, (calculated, given)
            %                 params.ambient.density,...         % two cases, (calculated, given)
            [   dtot_res, ...           % total damping coefficient
                ~, ...	% thermal damping coefficient
                ~, ...	% viscous damping coefficient
                ~] = ...          % polytropic exponent
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
            
            dtot_res = 2*dtot_res/(2*pi*fr); % total damping coefficient, M. Hall (1989)
            Delta_K_f_h = 1/(Para.rhow*pi*(ff^2))* r/((fr/ff)^2-1+1i*dtot_res);
            %         Delta_K_new = 1/(Para.rhow*pi*f^2)* r/((f0(indR0)/f)^2-1+1i*dtot_res(indR0))
            Delta_K = Delta_K + Delta_K_f_h; % compressibility of bubbles, M. Hall (1989)
            
        else
            continue
        end
        %         subplot(2,3,count),
        %     plot(xc*r*scale+x,yc*r*scale+z),hold on,axis equal,axis([-1 1 0 3])
    end
    % drawnow
    % % end
    % indR0;
    % figure(3)
    % R0 = sort(R0);
    % min(R0)
    % max(R0)
    % radii=[0.1e-4:0.1e-4:10e-3];
    % N = histcounts(R0,radii);
    % plot(N)
    % ylabel('bubble number')
    % xlabel('radii, mm')
    % title('radii, 1m')
    %
    % figure(4)
    % length(Fr)
    % % Fr = sort(Fr);
    % plot(Fr)
    % min(Fr)
    % max(Fr)
    % xlabel('bubble index')
    % ylabel('frequency, Hz')
    % title('Resonce frequency, 1m')
    %
    % figure(5)
    % fii=[1e3:1e3:1e5];
    % Nf = histcounts(Fr,fii);
    % plot(Nf)
    % xlabel('bubble number')
    % ylabel('frequency, kHz')
    % title('Resonce frequency, 1m')
    
    % for k=1:length(Tvalues)
    %     subplot(2,3,k),title([num2str(Tvalues(k)),' seconds']);
    % end
    
    
    waterVol = 1 * 0.8 * 0.9;
    Delta_K = Delta_K *1 /waterVol;
    Vg = Vg * 1/1;
    %     Vg;

    K_0 = 1/(Para.rhow*(Para.cw^2)); % compressibility of water

    Vf = Vg/waterVol; % Void fraction of the water column occupied by bubbles
%     Vf = 0;
    K = (1-Vf)*K_0 + Delta_K;  % compressibility of water with bubbles
%     Vf = 0;
%     K = K_0;
    params.effective.density = (1-Vf)*Para.rhow + Vf*params.gas.density; % effective density
    
    C_eff_m2 = K*params.effective.density; % Wood equation
    
%     C_eff_m2 = 1/(Para.cw^2) + Para.rhow*Delta_K; % 
    
    c_eff(countf) = sqrt(1/C_eff_m2); 
     
    c_b(countf) = real(c_eff(countf)); % Sound speed
    
    alpha_b(countf) = 20*2*pi*ff/log(10)*imag(1/c_eff(countf));  % Attenuation
    
    
end

% size(c_b)
% size(f)

figure(6)
plot(f*1e-3,c_b)
xlabel('frequency, kHz')
ylabel('sound speed, m/s')
% title(['Sound Speed, height = ',num2str(hh),'m'])
x0=200;
y0=50;
width=400;
height=240;
set(gcf,'units','points','position',[x0,y0,width,height])
% saveas(gcf,['SoundSpeedch4rate5Height',num2str(hh),'mnew.png'])
saveas(gcf, [tag, 'soundspeed.fig'])

figure(7)
plot(f*1e-3,abs(alpha_b))
xlabel('frequency, kHz')
ylabel('Attenuation, dB/m')
% title(['Attenuation, height = ',num2str(hh),'m'])
x0=200;
y0=50;
width=400;
height=240;
set(gcf,'units','points','position',[x0,y0,width,height])
% saveas(gcf,['Attenuationch4rate5Height',num2str(hh),'mnew.png'])
saveas(gcf, [tag, 'attenuation.fig'])

