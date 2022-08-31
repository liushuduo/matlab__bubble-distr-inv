
Nb=length(Bubbles);

scale=12;

Nbubs=length(Bubbles);

theta=linspace(0,2*pi,100);
xc=sin(theta);
yc=cos(theta);

% Tvalues=[100,200,300,400,500,600];
%Tvalues=20;

for k=1:Nbubs
    if (isempty(Bubbles(k).tstop)), Bubbles(k).tstop=inf; end
end

close
% for count=1:length(Tvalues)
%     T=Tvalues(count);

figure(6)
% filename = 'SoundSpeedChangeHeight.gif';
filename = 'BubbleCH4plume2Drate0p1dt1stime500swide.gif';
Tvalues = 0:1:500;
for count=1:length(Tvalues)
    disp(count);
    T=Tvalues(count);
%     R = [];
%     X = [];
%     Z = [];
    
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
%         R = [R r];
%         X = [X x];
%         Z = [Z z];
        
        %        subplot(1,1,count),plot(xc*r*scale+x*0.25,yc*r*scale+z),hold on,axis equal,axis([-1 1 0 3])
        
        %         subplot(2,3,count),
                plot(xc*r*scale+x,yc*r*scale+z),
                hold on,
                %         axis equal,
                axis([-3 3 0 120])
                xlabel('range, m')
                ylabel('height, m')
        %         ylim([0 120])
                title(['t=',num2str(T),'s'])
    end
%     plot(xc*R*scale+X,yc*R*scale+Z),
%     %         hold on,
%     %         axis equal,
%     axis([-2 2 0 120])
%     xlabel('range, m')
%     ylabel('height, m')
%     %         ylim([0 120])
%     title(['t=',num2str(T),'s'])
    
    drawnow
    x0=500;
    y0=10;
    width=200;
    height=700;
    set(gcf,'units','points','position',[x0,y0,width,height])
    set(gcf,'Renderer','Zbuffer');
    hold off
    %     scrnsz = get(0, 'ScreenSize');
    %     f = figure('Position', [1 1 100 100], 'MenuBar', 'none');
    frame = getframe(gcf);
    
    %     close(f);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if count == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
    end
    
end

% for k=1:length(Tvalues)
%     subplot(2,3,k),title([num2str(Tvalues(k)),' seconds']);
% end
