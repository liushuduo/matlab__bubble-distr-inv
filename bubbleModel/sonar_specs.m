function [Sonar]=sonar_specs(type)

% Function to define the structure Sonar containing the details of the specs for the sonars

switch(lower(type))
    
    case {'alds'}
        
        Sonar.type='alds';
        Sonar.range_res= 0.0375;
        Sonar.horz_res= 3.5*pi/180;
        Sonar.vert_res= 11*pi/180;
        Sonar.MRA_vert=0;
        Sonar.f=70000;

    case {'solstice'}
        
        Sonar.type='solstice';
        Sonar.range_res=0.0375;
        Sonar.horz_res= 0.15*pi/180;
        Sonar.vert_res=15*pi/180;
        Sonar.MRA_vert=-7.5*pi/180;
        Sonar.f=750000;
        
    otherwise
        error('Unknown sonar type')
end
