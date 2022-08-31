function [TS,sig]=plume_target_strength_struct(Bubbles,Time,dt,depth,f,depth_range)

if (exist('depth_range','var')~=1 || isempty(depth_range)), depth_range=[0 inf]; end

if numel(depth_range)~=2, error('depth_range must be a 2 element vector'), end

parameters_struct

Sig=0;
Nbubs=length(Bubbles);

sig=zeros(Nbubs,1);
for k=1:Nbubs
    
    if (isempty(Bubbles(k).tstop)),Bubbles(k).tstop=inf; end
    
    if ( Time>Bubbles(k).tstart && Time<Bubbles(k).tstop)
        
        m=round((Time-Bubbles(k).tstart)/dt)+1;
        if (m<=0), m=1;end
%        if (m>length(Bubbles(k).r)),keyboard,m=length(Bubbles(k).r); end
        %        disp([Bubbles(k).tstart Bubbles(k).tstop m])
        a=Bubbles(k).r(m);
        z=Bubbles(k).z(m);
        Pb=(Para.Patm+Para.rhow*Para.g*(depth-z)+2*Para.gamma_st/a);                                  % Equation (4) Leifer and Patro
        rhoCO2=Pb*Para.MCO2/(Para.R*(273+Para.T));
        
        if (z>depth_range(1) && z<depth_range(2))
            sig(k)  = cs_bubble_modal(a,Para.cCO2,Para.cw,rhoCO2,Para.rhow,f,pi);
        end
    end
    
end

TS=10*log10(sum(sig));
