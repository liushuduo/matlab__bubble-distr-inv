function [phi,phi1,phi2,phi3,Area]=leifer_culling_pop_minor(r)

% Seep bubble population from I. Leifer, D. Culling, Geo-Marine Letters, 30, 339-353, 2010, Figure 5a - minor vent
%
% INPUTS
%   r    Vector of radii, in m (paper quotes values in microns - conversion is performed in the function)
%
%  OUTPUTS
%  All populations are expressed as pdf's, i.e. they are normalised, and per second
%   phi                 Total population, which is a Gaussian mixture
%   phi1, phi2, phi3    Components of the Gaussian mixture: phi=phi1+phi2+phi3
% 

r=r*1e6;                     % Conversion from m to microns

phi1=0.01*exp( -((r-1878)/344).^2 );
phi2=0.0019*exp( -((r-2979)/436).^2 );
phi3=0.0004*exp( -((r-4149)/596).^2 );

phi=(phi1+phi2+phi3); 
Area=7.9881e-006;
phi=phi/Area;           % Normalise so area under curve is unity - conversion to pdf (7.9881e-006 computed off-line)   