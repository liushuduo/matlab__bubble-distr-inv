function [r,z,vb,t]=buoyant_bubble_v2(rstart,depth,dt)

% [r,z,vb,flux,t]=buoyant_bubble_v2(rstart,depth,dt)
%
% Implementation of bubble dynamics based on Leifer and Patro, Continental, Shelf Research, 2002, 2409-2428
% Uses a simple Euler approximation to integrate the differential equations
%
% INPUTS
% rstart - initial bubble radius [m]
% depth  - depth of bubble release [m]
% dt     - step size [s]
%
% OUTPUTS
% r      - vector of bubble radii [m]
% z      - vector of depths (above sea floor) [m]
% vb     - vector of bubble rise velocities [m/s]
% flux   - vector of molar flux [mol / s]
% t      - vector of time points [s]

parameters_struct

clean_flag=0;   % Dirty bubble
sea_flag=1;     % in seawater

count=1;

r(count)=rstart;
t(count)=0;
z(count)=0;
vb(count)=bubble_velocity_Fan_struct(r(count),clean_flag,sea_flag,Para);

while ( r(count)>1e-6 || z(count) > depth )  % Run the simulation until the bubble shrinks to 1 um

    [r(count+1),z(count+1),vb(count+1)]=bubble_updater_struct(r(count),z(count),vb(count),depth,dt,clean_flag,sea_flag,Para);
    
    t(count+1)=t(count)+dt;
    
    count=count+1;
end
