function [TS,sig]=plume_target_strength_shadow(Bubbles,Para,Sonar,Time,dt,depth,depth_range)

if (exist('depth_range','var')~=1 || isempty(depth_range)), depth_range=[0 inf]; end

if numel(depth_range)~=2, error('depth_range must be a 2 element vector'), end

Nbubs=length(Bubbles);

count=1;
for k=1:Nbubs
    
    if (isempty(Bubbles(k).tstop)),Bubbles(k).tstop=inf; end
    
    if ( Time>Bubbles(k).tstart && Time<Bubbles(k).tstop)       % Does bubble exist at the time of the TS calc?
        
        m=round((Time-Bubbles(k).tstart)/dt)+1;                 % Compute the time index corresponding to that time
        if (m<=0), m=1;end
        
        % Extracting info from Bubbles structure relating to only the bubbles at t=T
        A(count)=Bubbles(k).r(m);
        X(count)=Bubbles(k).x(m);
        Y(count)=Bubbles(k).y(m);
        Z(count)=Bubbles(k).z(m);
        count=count+1;
    end
end

Nbubs=length(A);        % Redefines Nbubs as the number of bubbles at the specified time

[Y,ind]=sort(Y,'descend');
A=A(ind);
X=X(ind);
Z=Z(ind);

ang=linspace(0,2*pi,360);
xc=cos(ang);
yc=sin(ang);

figure(1)
for k=1:Nbubs
    plot( xc*A(k)+X(k),yc*A(k)+Z(k)),hold on
end
axis equal

k=1;
while k<=length(A)     
    
    L=length(A);
    result=zeros(L,1);
    for n=k+1:L
        result(n)=shadow([X(k) Z(k)],A(k),[X(n) Z(n)],A(n));
    end
    
    ind=find(result);
    
    if ~isempty(ind)       
        A(ind)=[];
        X(ind)=[];
        Y(ind)=[];
        Z(ind)=[];
    end
    k=k+1;
end 

Nbubs=length(A);        % Redefines Nbubs as the number of bubbles after shadowed bubbles removed
figure(1)
for k=1:Nbubs
    plot( xc*A(k)+X(k),yc*A(k)+Z(k),'r'),hold on
end
axis equal

sig=zeros(Nbubs,1);
for k=1:Nbubs
    
    if (isempty(Bubbles(k).tstop)),Bubbles(k).tstop=inf; end
    
    if ( Time>Bubbles(k).tstart && Time<Bubbles(k).tstop)
        
        m=round((Time-Bubbles(k).tstart)/dt)+1;
        if (m<=0), m=1;end
%        if (m>length(Bubbles(k).r)),keyboard,m=length(Bubbles(k).r); end
        %        disp([Bubbles(k).tstart Bubbles(k).tstop m])
        a=A(k);
        z=Z(k);
        Pb=(Para.Patm+Para.rhow*Para.g*(depth-z)+2*Para.gamma_st/a);                                  % Equation (4) Leifer and Patro
        rhoCO2=Pb*Para.MCO2/(Para.R*(273+Para.T));
        
        if (z>depth_range(1) && z<depth_range(2))
            sig(k)  = cs_bubble_modal(a,Para.cCO2,Para.cw,rhoCO2,Para.rhow,Sonar.f,pi);
        end
    end
    
end

TS=10*log10(sum(sig));
