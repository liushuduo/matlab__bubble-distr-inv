function vb=bubble_velocity_Fan_struct(r,clean_flag,sea_flag,Para)

% function vb=bubble_velocity_Fan(r,clean_flag,sea_flag)
%
% INPUTS
% r          - Vector of bubble radii [m]
% clean_flag - 0,1 : 0 indicates 'dirty' bubble, 1 indicates clean bubble [Default: 1]
% sea_flag   - 0,1 : 0 indicates fresh water, 1 indicates sea water [Default: 1]
%
% OUTPUTS
% vb         - rise velocity [m / s]

% Based on "Bubble wake dynamics in liquids and liquid-solid suspensions", Book by Fan et al., as cited in
% Leifer & Patro Continental Shelf Research, 2002, 2409-2428, as (17) given
% there as equations 2-11 in Fan.

% Next 2 blocks commented out for speed
%if (exist('clean_flag','var')~=1 || isempty(clean_flag)), clean_flag=1; end
%if (exist('sea_flag','var')~=1 || isempty(sea_flag)), sea_flag=1; end

%if ( ~any(clean_flag==[0,1]) ), error('clean_flag must be 0 or 1'); end
%if ( ~any(sea_flag==[0,1]) ), error('sea_flag must be 0 or 1'); end

if (clean_flag)
    d=1.6;
else
    d=0.8;
end
if (sea_flag)
    c=1.4;
else
    c=1.2;
end

M=Para.g*Para.sigma_B^4/(Para.rhow*Para.gamma_st);          % Morton number
vb= ( (Para.rhow*Para.g*r.^2*M^0.038/(3.68*Para.sigma_B)).^(-d) + (c*Para.gamma_st./(Para.rhow*r) + Para.g*r).^(-d/2)  ).^(-1/d);
vb=real(vb);        % Occassional becomes complex (with small imaginary part) for small bubbles
