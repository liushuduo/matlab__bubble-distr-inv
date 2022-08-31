function [c,a,K,den,vp,vis,st]=waterprops(T,S,d,f)

% Returns  physical properties of water  for frequency f (Hz)
% at a given temperature T(oC), salinity S(ppt) and  depth d(m)
% nb f only used to calculate attenuation only but is still needed to run
% program
% outputs include:
% velocity c(m/s) 
% attenuation coefficient A(dB/m)
% bulk modulus K(Pa)
% density den(kg/m3)
% Vapour pressure vp(Pa) 
% Viscosity  vis(Pa.s) 
% surface tension st(N/m)

% this program calls subroutine cheb.m

% main reference is Siedler, P., Properties of seawater (pp 237-259 in
% "Numerical data and functional relationships in sccience and technology, ed. Sundermann, 1986)

% velocity and attenuation written by S. Arnott on 22/7/02,  
% other components added by Gary Robb on 20/09/05

%Calculation of velocity(m/s) from Lovett, JASA 63: 1713-1718.
% valid range T=-3 to 42 oC; S=0 to 41 ppt.
pc=1e-4*(1.01325e5+9.80665*d*1000);    % pressure in dbar required to get c
C3=1402.394;
CT=(5.01132*T)-(5.513036e-2*T^2)+(2.221008e-4*T^3);
CS=1.332947*S;
CP=(1.605336e-2*pc)+(2.12448e-7*pc^2);
CTSP=(-1.266383e-2*T*S)+(9.543664e-5*T^2*S)-(1.052396e-8*T*pc^2)+(2.183988e-13*T*pc^3)-(2.253828e-13*S*pc^3)+(2.062107e-8*T*S^2*pc);
c=C3+CT+CS+CP+CTSP;
clear pc C3 Ct CS CP CTSP

% calculation of bulk modulus of water (Pa) (Siedler, 1986 section 3.1.3)
% valid range: 0=<S=<42 (ppt); -2=<T=<40 (oC); 0=<pK=<1e4 (dbar)
pK=1e-4*(1.01325e5+9.80665*d*1000);   % pressure in dbar required to get K
K_pure=1.965221e5+1484.206*T-23.27105*T^2+1.360477e-1*T^3-5.155288e-4*T^4;
K_surface=K_pure+(546.746-6.03459*T+1.09987e-1*T^2-6.167e-4*T^3)*S+(7.944e-1+1.6483e-1*T-5.3009e-3*T^2)*S^1.5;
Aw=3.239908+1.43713e-3*T+1.16092e-4*T^2-5.77905e-7*T^3;
A=Aw+(2.2838e-3-1.0981e-5*T-1.6078e-6*T^2)*S+1.91075e-4*S^1.5;
Bw=8.50935e-6-6.12293e-7*T+5.2787e-9*T^2;
B=Bw+(-9.9348e-8+2.0816e-9*T+9.1697e-11*T^2)*S;
K_dbar=K_surface+A*pK+B*pK^2;
K=K_dbar*1e4;
clear pK k_pure K_surface Aw A Bw B

% calculation of density of water (kg/m3) (Siedler, 1986 section 3.1.3)
% valid range: 0=<S=<42 (ppt); -2=<T=<40 (oC); 0=<pK=<1e4 (dbar)
pd=1e-4*(1.01325e5+9.80665*d*1000);    % pressure in dbar required to get den
den_pure=999.842594+6.793952e-2*T-9.095290e-3*T^2+1.001685e-4*T^3-1.120083e-6*T^4+6.536332e-9*T^5;
den_surface=den_pure+(8.24493e-1-4.0899e-3*T+7.6438e-5*T^2-8.2467e-7*T^3+5.3875e-9*T^4)*S+(-5.72446e-3+1.0227e-4*T-1.6546e-6*T^2)*S^1.5+4.8314e-4*S^2;
den=den_surface/(1+(pd/K_dbar));
clear pd den_pure den_surface

%Calculation of attenuation (dB/m) from Schulkin and Marsh, JASA 34: 864-865.
% valid range =?
pa=1+9.80665*d*den/101325; % pressure for att calc in atmospheres
A1=2.34e-6;  
fT=21.9*10^(6-1520/(T+273));  
f1=f*1e-3;  %frequency (kHz)
B=3.38e-6;
a=8.686*(((S*A1*fT*f1^2)/(fT^2+f1^2))+((B*f1^2)/fT))*(1-(6.54e-4*pa));
clear pa A1 fT f1 B1

% Calculation of vapour pressure (siedler, 1986, section 3.1.8)
% ***omits any presure dependence***
% valid ranges: 0=<S=<40 (ppt); 0=<T=<40 (oC)
a_matrix=[1430.6181,-18.2465,7.6875,-0.0328,0.2728,0.1371,0.0629,0.0261,0.0200,0.0117,0.0067]; 
Tvp=T+273;  % temp in K
x=(2*Tvp-(648+273))/(648-273);
[chebx]=cheb(11,x);
logvp=(2794.0144/2+sum(a_matrix.*chebx(2:12)))/Tvp;
vp_pure=10^(3+logvp); % in Pa

A1=(-3.7433e-3)-(1.6537e-4*T)-(1.9667e-6*T^2)-(2.4350e-7*T^3);
B1=(5.2556e-4)-(7.7266e-6*T^3);
C1=(-4.9535e-5);
delta_vp=(A1*S)+(B1*S^1.5)+(C1*S^2);
vp=vp_pure+delta_vp;
clear A B C logvp chebx x Pvp a_matrix

%Calculation of viscosity vis(Pa.s) (siedler, 1986, section 3.1.12)
% valid range: 0=<S=<40 (ppt); 5=<T=<25 (oC) for surface water case
% pressure bit omitted as limited range of T and S and makes little difference
vis_pure=1.002e-3*10^((1.1709*(20-T)-1.827e-3*(T-20)^2)/(T+89.93));
cl=den*S/1806.55;
A=5.185e-5*T+1.0675e-4;
B=3.3e-5*T+2.591e-3;
vis=vis_pure*(1+A*cl^0.5+B*cl);

% Calculation of surface tension st (N/m) for distilled water
% (from Siedler, 1986, section 3.1.13)
% valid range: 0=<T=<60 (oC)
st=7.562e-2-1.3928e-4*T-3.063e-7*T^2;