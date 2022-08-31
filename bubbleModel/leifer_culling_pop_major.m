function [phi,phi1,phi2,phi3,phi4,Area]=leifer_culling_pop_major(r)

% Seep bubble population from I. Leifer, D. Culling, Geo-Marine Letters, 30, 339-353, 2010, Figure 5c - major vent
%
% INPUTS
%   r    Vector of radii, in m (paper quotes values in microns - conversion is performed in the function)
%
%  OUTPUTS
%  All populations are normalised, i.e. pdfs and expressed per second
%   phi                     Total population
%   phi1, phi2, phi3,phi4   Components of the Gaussian mixture: phi=phi1+phi2+phi3
% 

r=r*1e6;            % Conversion from m to microns

A1=(0.05*1000^2.22);
A2=(0.02*2000^0.21);

% Different components of the model for all r
phi1=(A1*r.^-2.22);
phi2=(A2*r.^-0.21);
phi3=0.0248*exp( -((r-3262)/329).^2 );
phi4=((0.006*5000^2.83)*r.^-2.83);

% Now compute the ranges over which each model applies

% Easy to explicitly compute the boundary between phi1 and phi2
rlim1=(A1/A2)^(1/2.01);
phi1=phi1.*(r<rlim1 & r>600);       % Plot in paper does not extend below 600 microns

% Have to approximate boundaries between phi2 and phi3 and phi3 and phi4
k=find(phi2<phi3,1,'first');
rlim2=r(k);
phi2=phi2.*(r>=rlim1 & r<rlim2);

k=find(phi4>phi3&(r>3262),1,'first');   % Look for crossing, but only after the peak in the Gaussian
rlim3=r(k);
phi4=phi4.*(r>=rlim3 & r<5000);         % Trend stops at 5000 microns

phi3=phi3.*(r>=rlim2 & r<rlim3);

phi=(phi1+phi2+phi3+phi4);
Area=1.06564e-004;
phi=phi/Area;         % Normalisation so result is a pdf, 1.06564e-004 is the integral's value across whole range
