function [TS,TS2D,Nb,Sb]=plume_target_strength_sonar_system(Bubbles,Para,Sonar,depth,Time,dt,range,plot_flag)

% Computes the target strength of the plume.  Computes value for resolution cells within the plume.
%
% [TS,sig]=plume_target_strength_sonar_system(Bubbles,Para,Sonar,depth,Time,dt,range)
%
% INPUTS:
% Bubbles    - Structure containing the info about the bubbles making up the plume
% Para       - Structure containing the physical parameters
% Sonar      - Structure containing the sonar system parameters
% depth      - Depth of the base of the plume [m]
% Time       - Time at which the plume's target strength is to be measured [s]
% dt         - Time step used in the generation of Bubbles [s]
% range      - Distance between the sonar and the plume [m]
% plot_flag  - Setting this plots the plume and resolution cells (default: 0, i.e. "off")
%
% OUTPUTS
% TS         - Target strength, 3 dimensional array containing the TS as a function of x,y and z [dB]
% TS2D       - Target strength summed vertically, matrix containing the TS in the x-y plane [dB]

if (exist('plot_flag','var')~=1 || isempty(plot_flag)), plot_flag=0; end

Nbubs=length(Bubbles);

Ehz=range*Sonar.horz_res/2;      % Ellipse's horizontal axis
Evt=range*Sonar.vert_res/2;      % Ellipse's vertical axis

Bub_count=0;

% Indentify the bubble which are present at the time T
count=1;
for k=1:Nbubs
    
    if (isempty(Bubbles(k).tstop)),Bubbles(k).tstop=inf; end
    
    if ( Time>Bubbles(k).tstart && Time<Bubbles(k).tstop)       % Does bubble exist at the time of the TS calc?
        
        m=round((Time-Bubbles(k).tstart)/dt)+1;                 % Compute the time index corresponding to that time
        if (m<=0), m=1;end
        
        % Extracting info from Bubbles structure relating to only the bubbles at t=T
        A(count)=Bubbles(k).r(m);
        X(count)=Bubbles(k).x(m);
        Y(count)=Bubbles(k).y(m);
        Z(count)=Bubbles(k).z(m);
        count=count+1;
    end
end

Nbubs=length(A);        % Redefines Nbubs as the number of bubbles at the specified time

if (plot_flag)
    close all
    figure(2)
    plot3(X,Y,Z,'.'),hold on
end

plume_height=max(Z);
plume_width=max(abs(X));
plume_depth=max(abs(Y));

% Computing the number of cells (range, azimuth and vertical) covering the plume
%Nheight=round(1+plume_height/(2*Evt));     % No overlap
Nheight=round(1+plume_height/Evt);          % 50% overlap

%N=ceil(plume_width/(2*Ehz));               % No overlap
N=ceil(plume_width/Ehz);                    % 50% overlap
Nwidth=2*N+1;

N=ceil(plume_depth/Sonar.range_res);        % No overlap (unshaded so overlap not used)
Ndepth=2*N+1;

[Xc,Zc,Yc]=cylinder(1,20);                  % Defines a cylinder on it side (hence odd ordering of the co-ordinates

Yc=(Yc-0.5)*Sonar.range_res;

Sb={};
disp([Nwidth Ndepth Nheight])
% Looping over the resolution cells in 3 for loops, indexed by p,q and r
for p=1:Nwidth
    %    Px= (p-(Nwidth+1)/2)*2*Ehz;            % No overlap
    Px= (p-(Nwidth+1)/2)*Ehz;              % 50% overlap
    
    for q=1:Nheight
        %        Pz=(q-1)*2*Evt;                    % No overlap
        Pz=(q-1)*Evt;                       % 50% overlap
        
        ang=linspace(0,2*pi,1000);
        Elx=Ehz*cos(ang);
        Elz=Evt*sin(ang);
        Ely=zeros(size(ang));
        if (plot_flag)
            surf(Xc*Ehz+Px,Yc,Zc*Evt+Pz,'FaceAlpha',0.5),colormap([1 0 0;1 0 0]),shading flat
            %plot3(Elx+Px,Ely,Elz+Pz,'r','LineWidth',2),hold on
        end
        
        for r=1:Ndepth
            
            Py= (r-(Ndepth+1)/2)*Sonar.range_res;
            % (Px,Py,Pz) should be the centre of the resolution cell, defined by the ellipse ( (x/Ehz)^2+(z/Evt)^2 = 1) and with depth Sonar.range_res
            
            % Is bubble in the resolution cell?
            Rad=((X-Px)/Ehz).^2 + ((Z-Pz)/Evt).^2;
            ind=find( Rad<1 & abs(Y-Py)<Sonar.range_res/2);
            
            Bub_count=Bub_count+length(ind);
            
            sig=zeros(Nbubs,1);
            
            % Compute the BSX for each bubble in the resolution cell
            for k=1:length(ind)
                % Bubble internal pressure
                Pb=(Para.Patm+Para.rhow*Para.g*(depth-Z(ind(k)))+2*Para.gamma_st/A(ind(k)));                                  % Equation (4) Leifer and Patro
                % Density of CO2 at pressure in bubble
                rhoCO2=Pb*Para.Mgas/(Para.R*(273+Para.T));
                
                % Compute the scattering cross-section
                sig(ind(k))  = cs_bubble_modal(A(ind(k)),Para.c_gas,Para.cw,rhoCO2,Para.rhow,Sonar.f,pi);
                if (isnan(sig(ind(k)))), keyboard, end

            end
            
            Sb{p,q,r}=sig;
            Nb(p,q,r)=length(ind);
            
            % Compute the overall target strength.  Note defintion of BSX used in Feuillade is such that TS=10 log10(BSX)
            S=sum(sig);
            if ( S == 0 )
                TS(p,q,r)=-200;
            else
                TS(p,q,r)=10*log10(sum(sig));
                
            end
        end
    end
end
ts=10.^(TS/10);
TS2D=10*log10(squeeze(sum(ts,2)));


