% % exportgraphics(fH, './figures/n_gt.png', 'Resolution', 600);
close all
% Get simulation result
f_sim = 100 : 100 : 9e3;
[c_b_sim, alpha_b_sim, c_eff_sim] = SoundspeedResponse(bubRadPop, f_sim, 1, params, Para);
measureFreq = [1 2 4 5 6 8] * 1e3;

c_eff_measure = [];
for f = measureFreq
    c_eff_measure(end+1) = c_eff_sim( f == f_sim );
end

% add noise 
sigma_n = 0.1;
c_eff_measure_n = c_eff_measure + sqrt(var(real(c_eff_sim))) * randn(1, length(c_eff_measure)) * sigma_n ...
                   + 1i * sqrt(var(imag(c_eff_sim))) * randn(1, length(c_eff_measure)) * sigma_n;
c_b_measure_n = real(c_eff_measure_n);
alpha_b_measure_n = - 2 * pi * 20 / log(10) * imag( 1 ./ c_eff_measure_n ) .* measureFreq;
alpha_b_measure_n( alpha_b_measure_n < 0 ) = 0;

% do interpolation
freqBin = 1e3 : 100 : 8e3;
c_b_inv = interp1(measureFreq, c_b_measure_n, freqBin, 'pchip');
alpha_b_inv = interp1(measureFreq, alpha_b_measure_n, freqBin, 'pchip');


% fig: simulation
palette = colororder;
fH = figure(Position=[1800 0 1000 800]);
tH = tiledlayout(1,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "inches";
tH.OuterPosition = [0 0 3.5 3.5/4*3];
aH = nexttile;
hold(aH,'on')

aH.YLim = [1.395, 1.56];
xlabel('frequency (kHz)')
ylabel('sound speed (km/s)')
sH = scatter(measureFreq/1e3, c_b_measure_n/1e3, 'LineWidth', 1);
sH.Marker = 'square';
sH.MarkerEdgeColor = '#00436f';
 sH.SizeData = 30;
lH = plot(freqBin/1e3, c_b_inv/1e3, 'LineStyle', '--', 'LineWidth', 1.5);
lH.Color = '#5bbeff';

yyaxis right

sH = scatter(measureFreq/1e3, alpha_b_measure_n, 'LineWidth', 1);
sH.MarkerEdgeColor = '#933811';
sH.Marker = 'diamond';
lH = plot(freqBin/1e3, alpha_b_inv, 'LineStyle', '--', 'LineWidth', 1.5);
lH.Color = '#ed8f67';
sH.SizeData = 30;

ylabel('attenuation (dB/m)')

aH.YAxis(1).Color = [0    0.4470    0.7410];
aH.YAxis(2).Color = [0.8500    0.3250    0.0980];
% aH.YAxis(1).Color = 'k';    aH.YAxis(2).Color = 'k';
grid(aH, 'on');
aH.FontName = 'Arial';
aH.FontSize = 7;
axis off

aH.YLim = aH.YLim + [-0.5, 0.5];
exportgraphics(fH, './figures/presentation_meas.png', 'Resolution', 600);

%% Get simulation result
clear; close all; clc
rng(6)
exportfig = true;

% Set up parameters
params = SetupParams('air');
Para = SetupEnvPara('air');

% All bubble radii for current simulation
load bubbledistr-experiment.mat
nBin = 21;
nBub = 1e3;
bubRadPop = randsample(bubRad, nBub, true, bubPdf);
% binEdge = linspace(min(bubRadPop), max(bubRadPop), nBin);
binEdge = linspace(0.2, 2.2, 21) * 1e-3;
n_gt = histcounts(bubRadPop, binEdge);


% fig: bubble_histogram
fH = figure(Position=[0, 0, 1000, 800]);
tH = tiledlayout(fH, 1, 1, 'TileSpacing','tight', Padding='tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 8 6];
aH = nexttile;
% palette = colorpalette('ieee_light');
% colororder(palette)
bH = bar(binEdge(1:end-1) * 1e3, n_gt);
bH.BarWidth =0.7;
bH.EdgeColor = "none";
% bH.FaceColor = palette{1};
aH.FontName = 'Arial';
aH.FontSize = 7;

% plot(bubRad*1e3,bubPdf/max(bubPdf)*max(n_gt)*1.1,'LineWidth', 1.5)
set(gca, "XTick", [])
xlabel("bubble radius")
aH.XLim = [0.2, 2];
axis off
exportgraphics(fH, './figures/presentation_n.png', 'Resolution', 600);



%% fig: bubble_histogram
fH = figure(Position=[0, 0, 1000, 800]);
tH = tiledlayout(fH, 1, 1, 'TileSpacing','tight', Padding='tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 8 6];
aH = nexttile;
% palette = colorpalette('ieee_light');
% colororder(palette)
plot(bubRad*1e3,bubPdf/max(bubPdf)*max(n_gt)*1.1,'LineWidth', 1.5)
% bH.FaceColor = palette{1};
aH.FontName = 'Arial';
aH.FontSize = 7;
axis off
exportgraphics(fH, './figures/presentation_Na.png', 'Resolution', 600);
%%
close all
% Get simulation result
f_sim = 100 : 100 : 9e3;
[c_b_sim, alpha_b_sim, c_eff_sim] = SoundspeedResponse(bubRadPop, f_sim, 1, params, Para);
measureFreq = [1 2 4 5 6 8] * 1e3;

c_eff_measure = [];
for f = measureFreq
    c_eff_measure(end+1) = c_eff_sim( f == f_sim );
end

% add noise 
sigma_n = 0.1;
c_eff_measure_n = c_eff_measure + sqrt(var(real(c_eff_sim))) * randn(1, length(c_eff_measure)) * sigma_n ...
                   + 1i * sqrt(var(imag(c_eff_sim))) * randn(1, length(c_eff_measure)) * sigma_n;
c_b_measure_n = real(c_eff_measure_n);
alpha_b_measure_n = - 2 * pi * 20 / log(10) * imag( 1 ./ c_eff_measure_n ) .* measureFreq;
alpha_b_measure_n( alpha_b_measure_n < 0 ) = 0;

% do interpolation
freqBin = 1e3 : 100 : 8e3;
c_b_inv = interp1(measureFreq, c_b_measure_n, freqBin, 'pchip');
alpha_b_inv = interp1(measureFreq, alpha_b_measure_n, freqBin, 'pchip');

fH = figure(Position=[1800 0 1000 800]);
tH = tiledlayout(1,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 8 6];
aH = nexttile;
hold(aH,'on')
lH = plot(f_sim/1e3, c_b_sim/1e3, 'LineWidth', 1);
lH.Color = [0    0.4470    0.7410];
aH.YLim = [1.395, 1.56];
xlabel('frequency (kHz)')
ylabel('sound speed (km/s)')


yyaxis right
lH = plot(f_sim/1e3, -alpha_b_sim, 'LineWidth', 1);
lH.Color = [0.8500    0.3250    0.0980];


ylabel('attenuation (dB/m)')

aH.YAxis(1).Color = [0    0.4470    0.7410];
aH.YAxis(2).Color = [0.8500    0.3250    0.0980];
% aH.YAxis(1).Color = 'k';    aH.YAxis(2).Color = 'k';
grid(aH, 'on');
aH.FontName = 'Arial';
aH.FontSize = 11;
axis off
exportgraphics(fH, './figures/presentation_th.png', 'Resolution', 600);


%% fig: simulation
palette = colororder;
fH = figure(Position=[1800 0 1000 800]);
tH = tiledlayout(1,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 12 8.3];
aH = nexttile;
hold(aH,'on')
lH = plot(f_sim/1e3, c_b_sim/1e3, 'LineWidth', 1);
lH.Color = [0    0.4470    0.7410];
aH.YLim = [1.395, 1.56];
xlabel('frequency (kHz)')
ylabel('sound speed (km/s)')
sH = scatter(measureFreq/1e3, c_b_measure_n/1e3, 'LineWidth', 1);
sH.Marker = 'square';
sH.MarkerEdgeColor = '#00436f';
 sH.SizeData = 30;
lH = plot(freqBin/1e3, c_b_inv/1e3, 'LineStyle', '--', 'LineWidth', 1);
lH.Color = '#5bbeff';

yyaxis right
lH = plot(f_sim/1e3, -alpha_b_sim, 'LineWidth', 1);
lH.Color = [0.8500    0.3250    0.0980];
sH = scatter(measureFreq/1e3, alpha_b_measure_n, 'LineWidth', 1);
sH.MarkerEdgeColor = '#933811';
sH.Marker = 'diamond';
lH = plot(freqBin/1e3, alpha_b_inv, 'LineStyle', '--', 'LineWidth', 1);
lH.Color = '#ed8f67';
sH.SizeData = 30;

ylabel('attenuation (dB/m)')

aH.YAxis(1).Color = [0    0.4470    0.7410];
aH.YAxis(2).Color = [0.8500    0.3250    0.0980];
% aH.YAxis(1).Color = 'k';    aH.YAxis(2).Color = 'k';
grid(aH, 'on');
aH.FontName = 'Arial';
aH.FontSize = 11;
axis tight
box on

legend('theoretical value', 'measurements', ...
       'interpolation', 'theoretical value', ...
       'measurements', 'interpolation', ...
        'NumColumns', 3, ...
       Location='northoutside', Orientation='horizontal')


aH.YLim = aH.YLim + [-0.5, 0.5];
exportgraphics(fH, './figures/presentation_fig_simulation.png', 'Resolution', 600);

%%
clear; close all; clc

% set bubble radii bin
nBubRadBin = 21;
bubRadBin = linspace(0.2, 2.2, nBubRadBin) * 1e-3;
bubRadList = bubRadBin(1 : end-1);

% get experiment data
f_exp = [1 2 4 5 6 8] * 1e3;
load experiment_data.mat


% recover complex sound speed
c_amb = mean(cExpMat(1, :));
iFlux = 4;                  % id 4 is 1 L/min
c_b_exp = cExpMat(iFlux, 1:length(f_exp));
alpha_b_exp = - (attenExpMat(iFlux, :) - attenExpMat(1, :));
alpha_b_exp( alpha_b_exp <= 0 ) = 0;
alpha_b_exp = alpha_b_exp(1:length(f_exp));
alpha_b_exp(end) = 5;

% do interpolation
freqBin = 1e3 : 100 : 8e3;
c_b_inv = interp1(f_exp, c_b_exp, freqBin, 'pchip');
alpha_b_inv = interp1(f_exp, alpha_b_exp, freqBin, 'pchip');

% recover complex sound speed
c_r = c_b_inv;
c_i = log(10) / 20 * alpha_b_inv .* c_r.^2 ./ (2*pi * freqBin);
c_eff_inv = c_r + 1i * c_i;

% fig:exp
fH = figure(Position=[1800 0 1000 800]);
tH = tiledlayout(1,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 12 8.3];
aH = nexttile;
hold(aH,'on')
xlabel('frequency (kHz)')
ylabel('sound speed (km/s)')
sH = scatter(f_exp/1e3, c_b_exp/1e3-0.13, 'LineWidth', 1);
sH.Marker = 'square';
sH.MarkerEdgeColor = '#00436f';
 sH.SizeData = 30;
lH = plot(freqBin/1e3, c_b_inv/1e3-0.13, 'LineStyle', '--', 'LineWidth', 1);
lH.Color = '#5bbeff';
aH.YLim = [1.44, 1.57]

yyaxis right
sH = scatter(f_exp/1e3, alpha_b_exp, 'LineWidth', 1);
sH.MarkerEdgeColor = '#933811';
sH.Marker = 'diamond';
lH = plot(freqBin/1e3, alpha_b_inv, 'LineStyle', '--', 'LineWidth', 1);
lH.Color = '#ed8f67';
sH.SizeData = 30;

ylabel('attenuation (dB/m)')

aH.YAxis(1).Color = [0    0.4470    0.7410];
aH.YAxis(2).Color = [0.8500    0.3250    0.0980];
% aH.YAxis(1).Color = 'k';    aH.YAxis(2).Color = 'k';
grid(aH, 'on');
aH.FontName = 'Arial';
aH.FontSize = 12;
axis tight

% legend('sound speed measurements','$\mathbf{c}_b$', ...
%     'attenuation measurements', '$\boldsymbol{\alpha}_b$', 'Interpreter', 'latex', ...
%     'NumColumns', 2, ...
%        Location='northoutside', Orientation='horizontal')

legend('sound speed measurements', 'interpolation', ...
       'attenuation measurements', 'interpolation', ...
        'NumColumns', 2, ...
       Location='northoutside', Orientation='horizontal')

box on
xlim([0, 9])
aH.YLim = aH.YLim + [-0.5, 0.5];
exportgraphics(fH, './figures/presentation_fig_exp.png', 'Resolution', 600);