%% Model Parameter defintion

% Basic parameter definitions
Duration=500;            % Length of simulation [s]
dt=1;                 % Time step [s]
depth=120;               % Depth of the base of the plume [m]
size_thresh=1e-6;       % Bubbles shrinking to this size are removed from the simulation [m]
sigma=0.03;              % Standard deviation of the random walk process to model horizontal migration of bubbles
Flux_rate=0.1;            % L / min at depth (not at the surface) later converted to m^3/s

% Setting flags to indicate clean bubble or dirty bubble and fresh or salt water
clean_flag=1;
sea_flag=1;

% Selecting the bubble size distribution model
model='chi';

% Generating the structure "Para" which contains the physical parameters
parameters_struct

% Bubble radii at which the pdf is evaluated so that samples can be drawn from it
radii=0.1e-3:0.1e-3:10e-3;

%% Computing derived parameters

Nits=Duration/dt+1;     % Number of time points

Flux_rate=Flux_rate/(1000*60);  % Convert to m^3/s

switch model
    case{'lcminor'}
        pdf=leifer_culling_pop_minor(radii);
    case{'lcmajor'}
        pdf=leifer_culling_pop_major(radii);
    case{'chi'}
        a=6;
        Rmean=2.6e-3;
        theta=(Rmean*gamma(a/2)/gamma(a/2+0.5)).^2;
        pdf=chi_pdf(radii,theta,a);
    case{'lognormal'}
        Rmean=2.6e-3;
        sigma0=7.65e-4;      % Mean and standard deviation set to approximately match chi distribution - which is a rather arbitrary choice
        pdf=lognormal_dist(radii,Rmean,sigma0);
    otherwise
        error('Unknown Model')
end

% Average bubble volume for the pdf
mean_vol=trapz(radii,4*pi*radii.^3.*pdf/3);
% Mean number of bubbles per second needed to generate required flux
Bubbles_per_sec=Flux_rate/mean_vol;

Nbub_it=Bubbles_per_sec*dt;             % Average number of bubbles starting in a iteration
% Nbub_it is unlikely to be an integer, so N1<N<N2 where N1 and N2 are integers.  For each time step N1 bubbles are
% generated with probability N-N1 (which should be a value in [0,1]) and N2 are generated with probability N2-N
% The effect is that the average number of bubbles generated per second is N (even if it is not an integer)
N1=floor(Nbub_it);
N2=ceil(Nbub_it);
prob=Nbub_it-N1;

% Time points for iterations
Times= 0:dt:Duration;

%% Preamble
% Nulling variables
index=[];
r=[];
vb=[];
z=[];
Bubbles=[];

% Setting up other variables
count=1;
% Imax=0;
h=waitbar(0,'Processing ..... ');
Ntot=0;

%% Looping over time instances
for t=Times
    waitbar(t/Times(end),h);
    if (rand>prob)
        N=N1;
    else
        N=N2;
    end
    Ntot=Ntot+N;
    
    % Generate N new bubbles from the pdf describing the bubble size distribution
    if N>0
        rnew=generate_rv(N,radii,pdf);
        znew=zeros(N,1);
        vbnew=bubble_velocity_Fan_struct(rnew,clean_flag,sea_flag,Para);
    else
        rnew=[];
    end
    
    % Update existing bubbles
    for k=index'
        [r,z,vb]=bubble_updater_struct(Bubbles(k).r(end),Bubbles(k).z(end),Bubbles(k).vb(end),depth,dt,clean_flag,sea_flag,Para);
        Bubbles(k).r(end+1) =r;
        Bubbles(k).vb(end+1)=vb;
        Bubbles(k).z(end+1) =z;
        steps=sigma*randn(2,1)*sqrt(dt)+[Para.c_slope*z;0]*dt;
        Bubbles(k).x(end+1) = Bubbles(k).x(end)+steps(1);
        Bubbles(k).y(end+1) = Bubbles(k).y(end)+steps(2);
        %        Bubbles(k).x(end+1) = Bubbles(k).x(end)+sigma*randn*dt+c_slope*z*dt;
        %        Bubbles(k).y(end+1) = Bubbles(k).y(end)+sigma*randn*dt;
    end
    
    % Remove bubbles which are too small
    m=1;
    remove_flag=zeros(size(index));
    for k=index'
        if (Bubbles(k).r(end)<size_thresh), remove_flag(m)=1; end
        m=m+1;
    end
    ind=find(remove_flag);
    for k=1:length(ind)
        Bubbles(index(ind(k))).tstop=t;
    end
    index(ind)=[];
    Imax = length(index);
    
    
    % Adding in the new bubbles
    for k=1:N
        Bubbles(Imax+k).r(1)=rnew(k);
        Bubbles(Imax+k).vb(1)=vbnew(k);
        Bubbles(Imax+k).z(1)=znew(k);
        Bubbles(Imax+k).x(1) =sqrt(rand)*cos(2*pi*rand);
        Bubbles(Imax+k).y(1) =sqrt(rand)*sin(2*pi*rand);
        Bubbles(Imax+k).tstart=t;
        Bubbles(Imax+k).tstop=inf;
    end
    index=[index;Imax+ (1:N)'];
%     Imax=Imax+N;
    count=count+1;
    
end
close(h)

save ('bubbleCH4plumerate0p1dt1stime500swide.mat','-v7.3');