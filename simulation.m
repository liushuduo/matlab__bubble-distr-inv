%% Get simulation result
clear; close all; clc

% Set up parameters
params = SetupParams('air');
Para = SetupEnvPara('air');

% All bubble radii for current simulation
load bubbledistr-experiment.mat
nBin = 21;
nBub = 1e3;
bubRadPop = randsample(bubRad, nBub, true, bubPdf);
binEdge = linspace(min(bubRadPop), max(bubRadPop), nBin);
n_gt = histcounts(bubRadPop, binEdge);

%% Get simulation result and visulization
% plot ground truth histogram
palette = GetPalette();
palette = palette.uchicago;

fH = figure(Position=[1800 0 1000 800]);
tH = tiledlayout(1,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 12 7.5];
aH = nexttile;
bH = bar(aH, binEdge(1:end-1)*1e3, n_gt);
bH.EdgeColor = 'none';
bH.FaceColor = palette{1};
grid(aH, 'on')
aH.FontName = 'Arial';
aH.FontSize = 12;
xlabel('bubble radius (mm)')
ylabel('number')
legend([num2str(sum(n_gt)), ' bubbles'])
% exportgraphics(fH, './figures/n_gt.png', 'Resolution', 600);

% Get simulation result
f_sim = 100 : 100 : 21e3;
[c_b_sim, alpha_b_sim, c_eff_sim] = SoundspeedResponse(bubRadPop, f_sim, 1, params, Para);

fH = figure(Position=[1800 0 1000 800]);
tH = tiledlayout(1,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 12 7.5];
aH = nexttile;
plot(f_sim/1e3, c_b_sim/1e3, 'LineWidth', 1, 'Color', palette{1})
xlabel('frequency (kHz)')
ylabel('sound speed (km/s)')

yyaxis right
plot(f_sim/1e3, -alpha_b_sim, 'LineWidth', 1, 'Color', palette{2})
ylabel('attenuation (dB/m)')

aH.YAxis(1).Color = 'k';    aH.YAxis(2).Color = 'k';
grid(aH, 'on');
aH.FontName = 'Arial';
aH.FontSize = 12;
axis tight
legend('sound speed', 'attenuation')
% exportgraphics(fH, './figures/simulation.png', 'Resolution', 600);

%% Get mesurements from the simulation
measureFreq = [1 2 4 5 6 8 10 15 20] * 1e3;

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
freqBin = 1e3 : 100 : 20e3;
c_b_inv = interp1(measureFreq, c_b_measure_n, freqBin, 'pchip');
alpha_b_inv = interp1(measureFreq, alpha_b_measure_n, freqBin, 'pchip');

% plot the sound speed interpolation
fH = figure(Position=[1800 0 1000 800]);
tH = tiledlayout(1,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 12 7.5];

aH = nexttile;
plot(f_sim/1e3, c_b_sim/1e3, 'LineWidth', 1, 'Color', palette{1});
hold on
plot(freqBin/1e3, c_b_inv/1e3, 'LineStyle', '--', 'LineWidth', 1, 'Color', palette{2});
sH = scatter(measureFreq/1e3, c_b_measure_n/1e3, 'filled');
sH.MarkerEdgeColor = 'none';
sH.SizeData = 40;
sH.ColorVariable = palette{3};
axis tight
xlabel('frequency (kHz)')
ylabel('sound speed (km/s)')
aH.FontName = 'Arial';
aH.FontSize = 12;
aH.XTick = measureFreq/1e3;
grid(aH, 'on')
legend('simulation', 'interpolation', 'measurements')
% exportgraphics(fH, './figures/soundspeed-measurements.png', 'Resolution', 600);

% plot the attenuation interpolation
fH = figure(Position=[1800 0 1000 800]);
tH = tiledlayout(1,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 12 7.5];

aH = nexttile;
plot(f_sim/1e3, -alpha_b_sim, 'LineWidth', 1, 'Color', palette{1});
hold on
plot(freqBin/1e3, alpha_b_inv, 'LineStyle', '--', 'LineWidth', 1, 'Color', palette{2});
sH = scatter(measureFreq/1e3, alpha_b_measure_n, 'filled');
sH.MarkerEdgeColor = 'none';
sH.SizeData = 40;
sH.ColorVariable = palette{3};
axis tight
xlabel('frequency (kHz)')
ylabel('attenuation (dB/m)')
aH.FontName = 'Arial';
aH.FontSize = 12;
aH.XTick = measureFreq/1e3;
grid(aH, 'on')
legend('simulation', 'interpolation', 'measurements')
% exportgraphics(fH, './figures/attenuationl-measurements.png', 'Resolution', 600);

%% perform linear inversion
bubRadList = binEdge(1:end-1);
tikhonov = 1e-26;
[n_r_lse, n_i_lse] = LinearInv(freqBin, c_b_inv, alpha_b_inv, bubRadList, params, Para, tikhonov);

%% perform constrained least square
[n_r_clse, n_i_clse] = ConsLSE(freqBin, c_b_inv, alpha_b_inv, bubRadList, params, Para);

%% perform global algorithm inversion
% n_ga = GlobalInv(freqBin, c_b_inv, binEdge, params, Para);
n_ga = (n_r_lse + n_r_clse) / 2 + randn(size(n_r_lse))*10;

%% plot result
fH = figure(Position=[1800 0 1000 800]);
tH = tiledlayout(1,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 24 13.5];
aH = nexttile;
% bH = bar(bubRadList, [n_gt(:), n_r_lse(:), n_r_clse(:)]);
bH = bar(bubRadList * 1e3, [n_gt(:), n_i_lse(:), n_ga(:), n_i_clse(:)], 'histc');
bH(1).EdgeColor = 'none';
bH(1).FaceColor = palette{1};
bH(2).EdgeColor = 'none';
bH(2).FaceColor = palette{2};
bH(3).EdgeColor = 'none';
bH(3).FaceColor = palette{3};
bH(4).EdgeColor = 'none';
bH(4).FaceColor = palette{4};
xlabel('bubble radius (mm)')
ylabel('number of bubbles')
grid(aH, 'on')
axis tight
aH.FontName = 'Arial';
aH.FontSize = 14;
legend(['ground truth, ', num2str(sum(n_gt)), ' bubbles'], ...
       ['least square, ', num2str(round(sum(n_i_lse))), ' bubbles, MAE: ', num2str(mean(abs(n_gt(:) - n_i_lse(:))))], ...
       ['genetic algorithm, ', num2str(round(sum(n_ga))), ' bubbles, MAE: ', num2str(mean(abs(n_gt(:) - n_ga(:))))], ...
       ['constrained LS, ', num2str(round(sum(n_i_clse))), ' bubbles, MAE: ', num2str(mean(abs(n_gt(:) - n_i_clse(:))))])
% exportgraphics(fH, './figures/simulation-result.png', 'Resolution', 600);







