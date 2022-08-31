
rad=(1:0.25:10) *1e-3;
lrad = length(rad);
depth=20;
dt=0.1;

count=1;
% R={};Z={};T={};
R = cell{lrad,1};
Z = cell{lrad,1};
T = cell{lrad,1};

Tmax=0;

for i = 1:length(rad)
    rstart=rad(i);
    [r,z,vb,t]=buoyant_bubble_v2(rstart,depth,dt);
    R{count}=r;
    Z{count}=z;
    T{count}=t;
    if (length(t)>Tmax), Tmax=length(t); end
    count=count+1;
end
count=count-1;

for k=1:count
    R{k}(Tmax)=0;
    Z{k}(Tmax)=0;
    T{k}=T{count};
end

inc=200e-3;
xcent=inc*(1:count);
theta=linspace(0,2*pi,100);
xc=sin(theta);
yc=cos(theta);
close all
h = zeros(count,1);
for k=1:count
    h(k)=plot(xc*R{k}(1)+xcent(k),yc*R{k}(1)+Z{k}(1));hold on,axis equal
end

scale=12;
for k=1:10:Tmax
    for m=1:count
        xd=xc*R{m}(k)*scale+xcent(m);
        yd=yc*R{m}(k)*scale+Z{m}(k);
        set(h(m),'Xdata',xd,'Ydata',yd)
        axis([-0.2 8 0 depth])
    end
    drawnow
end