function [] = GeneratePlume(fluxRate, model, tag, depth, duration, bubGenAreaFunc, varargin)
%GENERATEPLUME simulate bubble plume according to bubble transfer mechanism.

% INPUTS:
% fluxRate  - bubble generation rate [ L/min ]
% duration  - simulated bubble generation time [ s ]
% model     - name of bubble size distribution model 
% tag       - saved .mat file tag
% depth     - bubble generated depth

% OPTIONAL:
% dt        - time step [s]
% sizeTh    - bubbles shinking to this size are removed from the simulation
% sigma     - standard deviation of the random walk process to model horizontal migration of bubbles
% cleanFlag - clean bubble (1) or dirty bubble (0)
% seaFlag   - sea water (1) or fresh water (0)
% customBubPDFFile  - custom bubble size pdf file, 

%   Reproduced by Liu Shuduo
%   Signal Space and Information System Lab, ISEE, ZJU
%   Email:      sliu35@zju.edu.cn

% check input
p = inputParser;
addParameter(p, 'dt', 1)        % default time step: 1 s
addParameter(p, 'sizeTh', 1e-6) % default size threshold: 0.001 mm
addParameter(p, 'sigma', 0.03)  % default standard deviation: 0.03
addParameter(p, 'cleanFlag', 1) % default: clean bubble
addParameter(p, 'seaFlag', 0)   % default: fresh water 
addParameter(p, 'radii', 0.1e-3:0.1e-3:10e-3)    % default: 0.1 mm - 10 mm
addParameter(p, 'bubPdf', 0)    % required parameter when using costum model
parse(p, varargin{:})

% Gnerating the structure "Para" which contains the physical parameters
parameters_struct

% struct bubble size distribution
switch model
    
    case{'lcminor'}
        radii = p.Results.radii;
        pdf = leifer_culling_pop_minor(radii);

    case{'lcmajor'}
        radii = p.Results.radii;
        pdf = leifer_culling_pop_major(radii);

    case{'chi'}
        radii = p.Results.radii;
        a=6;
        Rmean=1e-3;
        theta=(Rmean*gamma(a/2)/gamma(a/2+0.5)).^2;
        pdf=chi_pdf(radii,theta,a);

    case{'lognormal'}
        radii = p.Results.radii;
        Rmean=2.6e-3;
        sigma0=7.65e-4;      % Mean and standard deviation set to approximately match chi distribution - which is a rather arbitrary choice
        pdf=lognormal_dist(radii,Rmean,sigma0);

    case{'custom'}
        
        if p.Results.bubPdf == 0
            error('bubPdf must be provided when using custom model!')
        end
        
        if length(p.Results.bubPdf) ~= length(p.Results.radii)
            error('radii vector must have same length with bubPdf vector!')
        end

        pdf = p.Results.bubPdf;
        radii = p.Results.radii;

    otherwise
        error('Unknown Model!')
end

% Computing derived parameters
timeInstList = 0 : p.Results.dt : duration;
nTimeInst = length(timeInstList);       % total number of time instances
fluxRate = fluxRate / (1000*60);        % convert to m^3/s

% Average bubble volume for the pdf
meanVol = trapz( radii, 4*pi/3 * radii.^3 .* pdf );

% Mean number of bubbles per second needed to generate required flux
nBubPerSec = fluxRate / meanVol;

% Mean number of bubbles generated per time instance
nBubPerTimeInst = nBubPerSec * p.Results.dt;
N1 = floor(nBubPerTimeInst);    N2 = ceil(nBubPerTimeInst);
prob = nBubPerTimeInst - N1;

% Nulling variables
livBubIndex = [];                         % living bubble index
Bubbles = [];                       % all bubbles status 
bubTimeTab = cell(nTimeInst, 1);

% setting up other variables
h = waitbar(0, 'Processing .....');
nBubTot = 0;           % total bubble number

disp('Start processing ...')
% Looping over time instances
for iTimeInst = 1 : nTimeInst

    % update waitbar
    t = timeInstList( iTimeInst );
    waitbar( t / timeInstList(end), h );
    
    % number of generated bubble number at current time instance
    if (rand > prob)
        N = N1;
    else 
        N = N2;
    end
    nBubTot = nBubTot + N;

    % generate N new bubbles from th epdf describing the bubble size distribution
    if N > 0
        rnew = generate_rv(N, radii, pdf);      % vector of generated bubble radii
        znew = zeros(N, 1);
        vbnew = bubble_velocity_Fan_struct(rnew, p.Results.cleanFlag, p.Results.seaFlag, Para);   % vector of bubble rise velocity
    else
        rnew = [];
    end
    
    % Update existing bubbles
    for k = livBubIndex'
        [r, z, vb] = bubble_updater_struct(Bubbles(k).r(end), Bubbles(k).z(end), Bubbles(k).vb(end), depth, p.Results.dt, p.Results.cleanFlag, p.Results.seaFlag, Para);
        Bubbles(k).r(end+1) = r;
        Bubbles(k).vb(end+1) = vb;
        Bubbles(k).z(end+1) = z;
        steps = p.Results.sigma * randn(2, 1) * sqrt(p.Results.dt) + [Para.c_slope*z; 0]*p.Results.dt;
        Bubbles(k).x(end+1) = Bubbles(k).x(end) + steps(1);
        Bubbles(k).y(end+1) = Bubbles(k).y(end) + steps(2);
    end

    % Remove bubbles which are too small or above water surface
    m = 1;
    removeFlagList = zeros(size(livBubIndex));
    for k = livBubIndex'
        if (Bubbles(k).r(end) < p.Results.sizeTh) || (Bubbles(k).z(end) >= depth)
            removeFlagList(m) = 1; 
        end
        m = m + 1;
    end

    % Remove bubbles
    ind = find(removeFlagList);
    for k = 1:length(ind)
        Bubbles(livBubIndex(ind(k))).tstop = t;
    end
    livBubIndex(ind) = [];

    Imax = length(Bubbles);
    % Adding in the new bubbles
    for k = 1:N
        % new bubble radius info
        Bubbles(Imax+k).r(1) = rnew(k);

        % new bubble rise velocity info 
        Bubbles(Imax+k).vb(1) = vbnew(k);

        % new bubble position info
        Bubbles(Imax+k).z(1) = znew(k);
        
        % generate randomly according to input function
        [Bubbles(Imax+k).x(1), Bubbles(Imax+k).y(1)] = bubGenAreaFunc();
        
        % new bubble time info
        Bubbles(Imax+k).tstart = t;
        Bubbles(Imax+k).tstop = inf;

    end
     
    livBubIndex = [livBubIndex; Imax + (1:N)'];

    bubTimeTab{ iTimeInst } = livBubIndex;
end

% save results
if ismac
    dataRootPath = '/Volumes/DATA/MATLAB/matlab__soundspeed-analysis';
elseif ispc 
    dataRootPath = 'S:/MATLAB/matlab__soundspeed-analysis';
end
fullFileName = fullfile(dataRootPath, 'bubble-model', [tag, '.mat']);
disp(['Saving ', fullFileName, ' ...']);
close(h);   clear h
save(fullFileName, '-v7.3');
disp('Finished!')
end