function [dtot, dth, dvis, drad] = bubble_damping_properetti(f,R0,cl,cg,Dl,Dg,Kl,Kg,gamma,rhol,rhog,sv)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Determine bubble damping coefficients given frequencies, radii 
% % and side parameters using Prosperetti 1977 theory
% % 
% % Inputs:
% % f is the driving frequency vector [Hz]
% % R0 is the bubble radius vector [m]
% % cl is the liquid sound speed [m/s]
% % cg is the gas sound speed [m/s]
% % Dl is the [length(R0) 1] surrounding liquid thermal diffusivity vector [m2/s]
% % Dg is the gas thermal diffusivity [m2/s]
% % Kl is the surrounding liquid thermal conductivity [W/(m*K)]
% % Kg is the gas thermal conductivity [W/(m*K)]
% % gamma is the gas ratio of specific heat []
% % rhol is the surrounding liquid density [kg/m3]
% % rhog is the [length(R0) 1] gas density vector [kg/m3]
% % sv is the ambient shear viscosity [Pa.s]
% % 
% % Outputs: 
% % dtot is the total damping coefficient []
% % dth is the thermal damping coefficient []
% % dvis is the viscous damping coefficient []
% % drad is the radiation damping coefficient []
% % kappa is the polytropic index []
% % 
% % Output dimensions:
% % drad: [length(f),length(R0)]
% % dth: [length(f),length(R0)]
% % dvis: [length(f),length(R0)]
% % dtot: [length(f),length(R0)]
% % kappa: [length(f),length(R0)]
% % 
% % dependancies: none
% %	source: Andrea Prosperetti|Thermal effects and damping mechanisms in 
% % the forced radial oscillations of gas bubbles in liquids|Journal of 
% % the Acoustical Society of America|vol. 61|No. 1|p17-27|January 1977
% % authorship: Benoit Berges
% % first created: 09/2011
% % last modified: 26/09/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = 2*pi*f;

% kappa = zeros(length(f),length(R0));
dth = zeros(length(f),length(R0));
dvis = zeros(length(f),length(R0));
drad = zeros(length(f),length(R0));
dtot = zeros(length(f),length(R0));

for Rin = 1:1:length(R0)
    for fin = 1:1:length(f)
        
        % calculating thermal diffusion length
        lth = sqrt(Dg(Rin)/(2*w(fin)));
        
        % calculating G1, have to use G1 approximation to make it work, why?
        G1 = (lth/ (cg/w(fin)))^2;
        
        % calculating G2
        G2 = w(fin)*R0(Rin)^2/Dg(Rin);
        
        % calculating G3
        G3 = w(fin)*R0(Rin)^2/Dl;
        
        % calculating G4
        G4 = Kl/Kg*(1 + (1 + 1i)*sqrt(G3/2));
        
        % calculating H1
        H1 = 1i + G1 + sqrt((1i - G1)^2 + 4*1i*G1/gamma);
        
        % calculating H2
        H2 = 1i + G1 - sqrt((1i - G1)^2 + 4*1i*G1/gamma);
        
        % calculating I1
        I1 = sqrt(gamma*G2/2*(1i - G1 + sqrt((1i - G1)^2 + 4*1i*G1/gamma)));
        
        % calculating I2
        I2 = sqrt(gamma*G2/2*(1i - G1 - sqrt((1i - G1)^2 + 4*1i*G1/gamma)));
        
        % calculating J1
        J1 = I1*coth(I1) - 1;
        
        % calculating J2
        J2 = I2*coth(I2) - 1;
        
        % calculating phi
        phi = (G4*(H2 - H1) + J2*H2 - J1*H1)/...
            (G4*(J2*H1 - J1*H2)- J1*J2*(H2 - H1));
        
        % calculating polytropic index
%         kappa(fin,Rin) = 1/3*gamma*G1*G2*real(phi);
        
        % calculating additional effective thermal viscosity
        mu_th = 1/4*w(fin)*rhog(Rin)*R0(Rin)^2*imag(phi);
        
        % calculating viscous damping
        dvis(fin,Rin) = 2*sv/(rhol*R0(Rin)^2);
        
        % calculating thermal damping
        dth(fin,Rin) = 2*mu_th/(rhol*R0(Rin)^2);
        
        % calculating acoustic damping
        drad(fin,Rin) = 1/2*w(fin)*(w(fin)*R0(Rin)/cl)/(1 + (w(fin)*R0(Rin)/cl)^2);
        
        % calculating total damping
        dtot(fin,Rin) = drad(fin,Rin) + dth(fin,Rin) + dvis(fin,Rin);
        
    end
end