
Nb=length(Bubbles);

scale=12;

Nbubs=length(Bubbles);

theta=linspace(0,2*pi,100);
xc=sin(theta);
yc=cos(theta);

% Tvalues=[100,200,300,400,500,600];
%Tvalues=20;
Tvalues = 1:4;

for k=1:Nbubs
    if (isempty(Bubbles(k).tstop)), Bubbles(k).tstop=inf; end
end

close
for count=1:length(Tvalues)
    T=Tvalues(count);
    for k=1:Nbubs
        
        if ( Bubbles(k).tstart<T && T<Bubbles(k).tstop )
            m=round((T-Bubbles(k).tstart)/dt)+1;
            m=max([m 1]);
        else
            continue
        end
        
        r=Bubbles(k).r(m);
        x=Bubbles(k).x(m);
        z=Bubbles(k).z(m);
%        subplot(1,1,count),plot(xc*r*scale+x*0.25,yc*r*scale+z),hold on,axis equal,axis([-1 1 0 3])

        subplot(2,3,count),plot(xc*r*scale+x,yc*r*scale+z),hold on,axis equal,axis([-2 2 0 5])
    end
    drawnow
end

for k=1:length(Tvalues)
    subplot(2,3,k),title([num2str(Tvalues(k)),' seconds']);
end
