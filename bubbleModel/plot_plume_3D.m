
Nb=length(Bubbles);

scale=12;

Nbubs=length(Bubbles);

[xs,ys,zs]=sphere(20);

% Tvalues=[1,2,5,10,15,20];
%Tvalues=20;
Tvalues = 1:5;

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
        y=Bubbles(k).y(m);
        z=Bubbles(k).z(m);
%        subplot(1,1,count),plot(xc*r*scale+x*0.25,yc*r*scale+z),hold on,axis equal,axis([-1 1 0 3])

        % Scaling and shifting sphere info
        xp=xs*r*scale+x;
        yp=ys*r*scale+y;
        zp=zs*r*scale+z;
        
        subplot(2,3,count),surf(xp,yp,zp,'FaceAlpha',0.5),colormap([0 0 1;0 0 1]),shading flat,hold on,axis equal,axis([-1 1 -1 1 0 3])
    end
    drawnow
end

for k=1:length(Tvalues)
    figure(1)
    subplot(2,3,k),title([num2str(Tvalues(k)),' seconds']);
end
