
% define_phys_parameters;
% load Bubbles120m0p01Lminxh.mat
clear all
close all

load S:\MATLAB\matlab__soundspeed-analysis\bubble-model\zhoushan_flux-rate-0.2_bubpdf-experiment.mat
% load bubbleCH4plumerate5dt10stime500swide.mat;

% Basic parameter definitions
% Duration=500;            % Length of simulation [s]
% dt=0.1;                 % Time step [s]
% depth=120;               % Depth of the base of the plume [m]
% size_thresh=1e-6;       % Bubbles shrinking to this size are removed from the simulation [m]
% sigma=0.03;              % Standard deviation of the random walk process to model horizontal migration of bubbles
% Flux_rate=0.01;            % L / min at depth (not at the surface) later converted to m^3/s

parameters_struct;

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

% Nsignals = 201;
% figure(1);
% filename = 'SoundSpeedChange3kHzTimerate0p2h2ra0p2.gif';
% Symbol = 200; % No. of symbol of the signal
% symbolTag = '200';


ff = 3000; % acoustic frequency
% f = 100:100:10000;
% for countf = 1:length(f)
%     countf
%     ff = f(countf);

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

%     Tvalues=[1,2,5,10,15,20];
% Tvalues=0:10:500;
Tvalues = 1:5;

for k=1:Nbubs
    if (isempty(Bubbles(k).tstop)), Bubbles(k).tstop=inf; end
end

close

h = 0:.5:5;
ra = -0.5:0.1:0.5;

Vf = zeros(length(Tvalues),length(h),length(ra));
C_eff = zeros(length(Tvalues),length(h),length(ra));
C_b = zeros(length(Tvalues),length(h),length(ra));
Alpha_b = zeros(length(Tvalues),length(h),length(ra));

for count=1:length(Tvalues)
    disp(count);
    T=Tvalues(count);
    
    %     T = 20;
    
    % h = 0.6;
    % ra = -0.1;
    %     vg = 0;
    
    
    %     Nob = zeros(length(h),length(ra));
    
    %     kr = [];
    vf = zeros(length(h),length(ra));
    c_eff = zeros(length(h),length(ra));
    c_b = zeros(length(h),length(ra));
    alpha_b = zeros(length(h),length(ra));
    
    for counth = 1:length(h)
        %         counth;
        hh = h(counth);
        for countra = 1:length(ra)
            rra = ra(countra);
            % for h = 0.1:0.1:3
            %             indR0 = 0;
            %             R0 = [];
            %             Fr = [];
            Delta_K = 0;
            % figure(100)
            vg = 0;
            %             X = [];
            %             Nozx = 0;
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
                %                 X(k) = x;
                %        subplot(1,1,count),plot(xc*r*scale+x*0.25,yc*r*scale+z),hold on,axis equal,axis([-1 1 0 3])
                if ( z <= hh && z > hh - 1 )
                    %                     z;
                    if x <= rra + 0.05 && x > rra - 0.05
                        x;
                        %                         Nozx = Nozx +1;
                        %                         indR0 = indR0 +1;
                        %                         R0(indR0) = r;
                        %                         kr = [kr k];
                        vg = vg + 4/3 * pi * r^3;  % Gas volume
                        %vg                         r;
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
                        %                         fr;
                        % calculating resonant frequencies using leighton model
                        %         fr = ...                         % resonant frequencies
                        %             bubble_resonance_leighton(  r,...
                        %             params.ambient.shear_viscosity,...
                        %             params.gas.ratio_spec_heat,...
                        %             params.ambient.vapor_pressure,...
                        %             params.ambient.surface_tension,...
                        %             params.ambient.pressure,...
                        %             params.ambient.density,...
                        %             1);
                        %                         fr;
                        %                         Fr = [Fr fr];
                        %         size(fr)
                        %         pause
                        %                         f0(indR0) = fr;
                        %                 params.ambient.sound_speed,...     % two cases, (given, calculated)
                        %                 params.gas.sound_speed,...         % two cases, (calculated, given)
                        %                 params.ambient.density,...         % two cases, (calculated, given)
                        [   dtot_res, ...           % total damping coefficient
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
                        
                        dtot_res = 2*dtot_res/(2*pi*fr); % M. Hall (1989)
                        Delta_K_f_h = 1/(Para.rhow*pi*(ff^2))* r/((fr/ff)^2-1+1i*dtot_res);
                        %         Delta_K_new = 1/(Para.rhow*pi*f^2)* r/((f0(indR0)/f)^2-1+1i*dtot_res(indR0))
                        Delta_K = Delta_K + Delta_K_f_h; % compressibility of bubbles, M. Hall (1989)
                    else
                        continue
                    end
                else
                    continue
                end
                %         subplot(2,3,count),
                %             plot(xc*r*scale+x,yc*r*scale+z),hold on,axis equal,axis([-1 1 0 3])
                
            end
            %             Nob(counth,countra) = Nozx;
            %         drawnow
            % end
            % indR0;
%             figure(3)
%             R0 = sort(R0);
%             min(R0)
%             max(R0)
%             radii=[0.1e-4:0.1e-4:10e-3];
%             N = histcounts(R0,radii);
%             plot(N,'.-')
%             ylabel('bubble number')
%             xlabel('radii, mm')
%             title(['radii, height = ',num2str(h),'m'])
            
%             figure(4)
%             length(fr)
%             % Fr = sort(Fr);
%             plot(fr)
%             min(fr)
%             max(fr)
%             xlabel('bubble index')
%             ylabel('frequency, Hz')
%             title(['Resonce frequency, height = ',num2str(h),'m'])
%             
%             figure(5)
%             fii=[1e3:1e3:1e5];
%             Nf = histcounts(Fr,fii);
%             plot(Nf)
%             xlabel('bubble number')
%             ylabel('frequency, kHz')
%             title(['Resonce frequency, height = ',num2str(h),'m'])
%             
%             for k=1:length(Tvalues)
%                 subplot(2,3,k),title([num2str(Tvalues(k)),' seconds']);
%             end
            Delta_K = Delta_K * 1^3/(1*0.1*4);
%             Vg;
            
            K_0 = 1/(Para.rhow*Para.cw^2);
%                     Vg
            vf(counth,countra) = vg/(1*0.1*4); % Void fraction of the water column
                    if hh == 0.6
                        disp(Vf)
                    end
            K = (1-vf(counth,countra))*K_0 + Delta_K;
            
            params.effective.density = (1-vf(counth,countra))*Para.rhow + vf(counth,countra)*params.gas.density;
            
            C_eff_m2 = K*params.effective.density; % Wood equation
%             
%                         C_eff_m2 = 1/(Para.cw^2) + Para.rhow*Delta_K; %
            
%             c_eff = sqrt(1/C_eff_m2);
            
%             c_b = real(c_eff);
            
%             alpha_b = 20*2*pi*ff/log(10)*imag(1/c_eff);
            c_eff(counth,countra) = sqrt(1/C_eff_m2);
            
            c_b(counth,countra) = real(c_eff(counth,countra)); % Sound speed
            
            alpha_b(counth,countra) = 20*2*pi*ff/log(10)*imag(1/c_eff(counth,countra));  % Attenuation
%             figure(6)
%             plot(X)
        end
    end
    
%     size(c_b)1
% %     size(f)
%     
%         figure(6)
%         surf(ra,h-2, c_b)
%         view(0, 90);
%         c = colorbar;
%         c.Label.String = 'Sound Speed, m/s';
%         % c.Limits = [1450 1550];
%         caxis([1460, 1540]);
%         colormap jet
%         shading interp
%         ylabel('height, m')
%         ylim([0 120])
%         xlabel('range, m')
%         zlabel('sound speed, m/s')
%         zlim([1460 1540])
%         title(['Sound Speed, T =',num2str(T),'s'])
%         %     saveas(gcf,['SoundSpeedHeightRange',num2str(T),'s.png'])
%         drawnow;
%         %     set(gcf,'outerposition',get(0,'screensize'));
%         x0=400;
%         y0=10;
%         width=250;
%         height=500;
%         set(gcf,'units','points','position',[x0,y0,width,height])
%         set(gcf,'Renderer','Zbuffer');
%         %     scrnsz = get(0, 'ScreenSize');
%         %     f = figure('Position', [1 1 100 100], 'MenuBar', 'none');
%         frame = getframe(gcf);
%     
%         %     close(f);
%         im = frame2im(frame);
%         [imind,cm] = rgb2ind(im,256);
%         if count == 1
%             imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.2);
%         else
%             imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
%         end
%     
%     
%         figure(7)
%         surf(ra,h-1,abs(alpha_b))
%         view(0, 90);
%         c = colorbar;
%         c.Label.String = 'Attenuation, dB/m';
%         caxis([0, 4]);
%         colormap jet
%         shading interp
%         ylabel('height, m')
%         ylim([0 1.0])
%         xlabel('range, m')
%         zlabel('Attenuation, dB/m')
%         zlim([0 4])
%         title(['Attenuation, T =',num2str(T),'s'])
%         saveas(gcf,['AttenuationHeightRange',num2str(T),'s.png'])
%         drawnow;
%         %     set(gcf,'outerposition',get(0,'screensize'));
%         x0=10;
%         y0=10;
%         width=500;
%         height=600;
%         set(gcf,'units','points','position',[x0,y0,width,height])
%         set(gcf,'Renderer','Zbuffer');
%         %     scrnsz = get(0, 'ScreenSize');
%         %     f = figure('Position', [1 1 100 100], 'MenuBar', 'none');
%         frame = getframe(gcf);
%     
%         %     close(f);
%         im = frame2im(frame);
%         [imind,cm] = rgb2ind(im,256);
%         if count == 1
%             imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.2);
%         else
%             imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
%         end
%     
%     figure(8)
%     surf(ra,h-0.1,Vf)
%     c = colorbar;
%     c.Label.String = 'Void Fraction';
%     colormap jet
%     shading interp
%     ylabel('height, m')
%     ylim([0 2.4])
%     xlabel('range, m')
%     zlabel('Void Fraction')
%     title(['Void Fraction'])
%     saveas(gcf,['VoidFractionHeightRange',num2str(ff),'Hz.png'])
%     
%     figure(9)
%     surf(ra,h-0.1,Nob)
%     c = colorbar;
%     c.Label.String = 'Bubble number';
%     colormap jet
%     shading interp
%     ylabel('height, m')
%     ylim([0 2.4])
%     xlabel('range, m')
%     zlabel('Bubble number')
%     title(['Bubble number'])
%     saveas(gcf,['Bubble numberHeightRange',num2str(ff),'Hz.png'])
    Vf(count,:,:) = vf;
    C_eff(count,:,:) = c_eff;
    C_b(count,:,:) = c_b;
    Alpha_b(count,:,:) = alpha_b;
    
end

save SoundSpeedch4wideChange3kHzTime500sdt10srate5h1ra0p1_.mat Vf C_eff C_b Alpha_b
