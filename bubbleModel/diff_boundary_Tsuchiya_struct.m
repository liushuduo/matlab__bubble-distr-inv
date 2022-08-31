function [kbub,kbub1,kbub2]=diff_boundary_Tsuchiya_struct(r,vb,Para,clean_flag)

% function kbub=diff_boundary_Tsuchiya(r,vb)
%
% INPUTS
% r          - Vector of bubble radii [m]
% vb         - Velocity of rising bubble [m/s]
% Para       - Structure containing physical parameters
% clean_flag - Flag indicating whether bubbles are clean
%
% OUTPUTS,
% kbub      - gas transfer rate [m / s]

% Based on Leifer & Patro Continental Shelf Research, 2002, 2409-2428 and Tsuchiya et al., Chemical Engineering Phsyics, Vol. 52, Nos 21-22, 4119-4126
% which is also following Clift and others

Re=2*r.*vb/Para.nu;

% Convert quantities from m to cm before application of (19)-(22)
gcm=Para.g*100;
Dfcm=100^2*Para.Df;
nucm=100^2*Para.nu;
rcm=100*r;
vbcm=100*vb;

if (clean_flag)     % Clean bubbles
    kbub1=sqrt( (2/pi) * (1 - 2.89./sqrt(Re)) .* (vb*Para.Df./r) );
    
    kbub1(imag(kbub1)~=0)=0;
    
    kbub2=sqrt(0.212*Dfcm.*vbcm./rcm)/100;
    
    kbub=max([kbub1(:) kbub2(:)],[],2);
%     if kbub1>kbub2
%         kbub=kbub1;
%     else
%         kbub=kbub2;
%     end
    
else                % Dirty bubbles
    zeta=10.^(0.5*(tanh(3.9*log10(2000*r/0.87))-1));    % Based on (11) in Tsuchiya et al., Chemical Engineering Phsyics, Vol. 52, Nos 21-22, 4119-4126
    % where de is a diameter measured in mm
    
    kbub1=sqrt( (2/pi) * (1 - 2.89./sqrt(Re)) .* (Para.Df*vb.*zeta./r) );
    kbub1(imag(kbub1)~=0)=0;
    
    kbub2=(0.42*gcm^0.3*(Dfcm/nucm)^(2/3)*nucm^0.4*(rcm).^-0.1)/100;  % (22) Leifer and Patro
    kbub1(imag(kbub2)~=0)=0;
    
    %kbub3=sqrt(0.212*Para.Df*vb./r)/100;       % (19) Leifer and Patro - not used
    %kbub3=(0.5*vbcm+Dfcm/rcm)/100;       % (21) Leifer and Patro - does not look reasonable
    
    kbub=max([kbub1(:) kbub2(:)],[],2);
%     if kbub1>kbub2
%         kbub=kbub1;
%     else
%         kbub=kbub2;
%     end
end