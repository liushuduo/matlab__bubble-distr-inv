function sig = cs_bubble_modal(a,c1,c,rho1,rho,f,theta)

% Compute the acoustic cross-section of a bubble, does not require a small ka assumption
%
%       sig = cs_bubble_modal(a,c1,c,rho1,rho,f,theta)
%
% INPUTS:
%   a       - Bubble radius [m]
%   c1      - Speed of sound in the gas in the bubble [m/s]
%   c       - Speed of sound in the fluid surrounding the bubble [m/s]
%   rho1    - Density of the gas in the bubble [kg/m^3]
%   rho     - Density of the fluid surrounding the bubble [kg/m^3]
%   f       - Vector of frequencies to evaluate the cross-section [Hz]
%   theta   - Scattering direction [Radians - pi = backscattered signal]
%
% OUTPUT:
%   sig    - Backscatter cross-section
%
% From Feuillade and Clay "Anderson (1950) Revisited" JASA, 106(2), 553-564, 1999

h = c1/c;
g = rho1/rho;

w = 2*pi*f;
k = w/c;
k1 = w/c1;

ka = k*a;
ka1 = k1*a;

Nmax=max(ka)+30;
% Nmax=max(ka)+100;
m=0:Nmax;

jm = sphbeslj(m,ka);
djm = sphbesldj(m,ka);
hm = sphhn(m,ka);
dhm=sphdhn(m,ka,1,hm);

jm1 = sphbeslj(m,ka1);
djm1 = sphbesldj(m,ka1);

Dm = (djm1.*jm-g*h*jm1.*djm)./(djm1.*hm-g*h*jm1.*dhm);      % Equation (9a) F&C

L = zeros(length(ka),1);

for indka = 1:1:length(ka)
    S = 0;
    for indm = 1:1:Nmax
        Pm = legendre(m(indm),cos(theta));
        Pm = Pm(1,1);
        S = S + Pm*(2*m(indm) + 1)*Dm(indka,indm);
    end
    L(indka) = 1i*a(indka)./ka(indka)*S;
end

sig = abs(L).^2;

end

% derivative of the spherical bessel function of the first kind
% function 	djn=sphbesl(n,x,flag,jn0)
%  flag = 0: compute jn0
%       = 1: do not compute jn0, with jn0 as input
%  djn(i)=i*jn(i)/x-jn(i+1)

function 	djn=sphbesldj(n,x,flag,jn0)

if nargin < 3
    flag=0;
end

indx=find(abs(x) < 1e-14);
x(indx)=1e-10*ones(size(indx));

x=x(:);
% m=length(x);
if ( flag == 0 )
    jn=sphbeslj(min(n):max(n)+1,x);
else
    jn_max=sphbeslj(max(n)+1,x);
    jn=[jn0 jn_max];
end

if (length(n) == 1)
    djn(:,1)=n*jn(:,1)./x-jn(:,2);
    return
else
    djn = zeros(size(jn,1),max(n)-min(n)+1);
    for i=1:max(n)-min(n)+1
        djn(:,i)=n(i)*jn(:,i)./x-jn(:,i+1);
    end
end

end

% derivative of the spherical bessel function of the second kind
% function 	dyn=sphbesldy(n,x,flag,jn0)
%  flag = 0: without yn0
%       = 1: do not compute yn0, with yn0 as input
%  dyn(i)=i*yn(i)/x-yn(i+1)

function 	dyn=sphbesldy(n,x,flag,yn0)

% flag = 0    : compute yn
% flag = 1    : do not compute yn

if nargin < 3
    flag=0;
end

x=x(:);
if ( flag == 0 )
    yn=sphbesly(min(n):max(n)+1,x);
else
    yn_max=sphbesly(max(n)+1,x);
    yn=[yn0 yn_max];
end

if (length(n) == 1)
    dyn=n*yn(:,1)./x-yn(:,2);
    return
else
    dyn = zeros(size(yn,1),max(n)-min(n)+1);
    for i=1:max(n)-min(n)+1
        dyn(:,i)=n(i)*yn(:,i)./x-yn(:,i+1);
    end
end
end

% spherical bessel function of the first kind
% function 	jn=sphbeslj(n,x)

function 	jn=sphbeslj(n,x)

indx=find(abs(x) < 1e-14);
x(indx)=1e-10*ones(size(indx));
x=x(:);
% size(n+1/2)
% size(x)
% Jm_1_2=besselj(n+1/2,x);
Jm_1_2 = zeros(length(x),length(n));
for k = 1:1:length(n)
    for l = 1:1:length(x)
        Jm_1_2(l,k)=besselj(n(k)+1/2,x(l));
    end
end

if (length(n) == 1)
    jn=sqrt(pi./(2*x)).*Jm_1_2;
    return
else
    jn = zeros(size(Jm_1_2,1),max(n+1)-min(n+1)+1);
    for i=1:max(n+1)-min(n+1)+1
        jn(:,i)=sqrt(pi./(2*x)).*Jm_1_2(:,i);
    end
end
end

% spherical bessel function of the second kind
% forward recurrence formula
% function      yn=sphbesly(n,x)

function      yn=sphbesly(n,x)

x=x(:);
Nmax=max(n+1);
Nmin=min(n+1);
yn = zeros(length(x),Nmax);
for j=1:Nmax
    if j==1
        yn(:,j)=-cos(x)./x;
    else
        if j==2
            yn(:,j)=-(cos(x)./x+sin(x))./x;
        else
            yn(:,j)=(2*j-3)*yn(:,j-1)./x-yn(:,j-2);
        end
    end
end
yn=yn(:,Nmin:Nmax);

end

% derivative of the spherical hankel function
% function dhn=sphdhn(n,x,flag,hn)
%  flag = 0: without yn0
%       = 1: do not compute yn0, with yn0 as input
%  dhn(i)=i*hn(i)/x-hn(i+1)

function dhn=sphdhn(n,x,flag,hn)

if nargin <= 3
    dhn=sphbesldj(n,x)+1i*sphbesldy(n,x);
elseif flag == 1
    dhn=sphbesldj(n,x,flag,real(hn))+1i*sphbesldy(n,x,flag,imag(hn));
end

end

% spherical hankel function
% function hn=sphhn(n,x)

function hn=sphhn(n,x)

hn=sphbeslj(n,x)+1i*sphbesly(n,x);
end
