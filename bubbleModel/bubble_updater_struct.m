function [r,z,vb]=bubble_updater_struct(r,z,vb,depth,dt,clean_flag,sea_flag,Para)

% Function to increment the bubble properties from time t to t+dt
%
% INPUTS
% r -       Bubble radius at time t [m]
% z -       Bubble height above sea floor at time t [m]
% vb -      Bubble's rise velocity at time t [m/s]
% depth -   Depth of the base of the bubble plume [m]
% dt -      Time step [s]
% clean_flag    Indicating whether bubble is clean 0 - dirty bubble, 1 - clean bubble
% sea_flag      Indicating fresh or sea water 0 - fresh water, 1 - sea water
% Para          Structure containing the physical parameters
%
% OUTPUTS
% r -       Bubble radius at time t+dt [m]
% z -       Bubble height above sea floor at time t+dt [m]
% vb -      Bubble's rise velocity at time t+dt [m/s]

% Based on Leifer and Patro Contential Shelf Research, 2002, Vol 22, 2409-2428

kbub=diff_boundary_Tsuchiya_struct(r,vb,Para,clean_flag);
Pb=(Para.Patm+Para.rhow*Para.g*(depth-z)+2*Para.gamma_st/r);                                  % Equation (4) Leifer and Patro
% Para.C-Pb/Para.H;
flux=kbub*4*pi*r^2*(Para.C-Pb/Para.H);
% flux = 0;
q=3*Para.R*(273+Para.T)/(4*pi*r^3);
% q=3*Para.R*(Para.T)/(4*pi*r^3);
dr=dt*(r*(q*flux-Para.rhow*Para.g*vb)/( 3*(Para.Patm+Para.rhow*Para.g*(depth-z))+4*Para.gamma_st/r ));
% dr = 0;
% Equation (8) Leifer and Patro
r=r+dr;
vbold=vb;
vb=bubble_velocity_Fan_struct(r,clean_flag,sea_flag,Para);
z=z+(vbold+vb)*dt/2;