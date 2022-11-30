clear; close all; clc

% set bubble radii bin
nBubRadBin = 21;
bubRadBin = linspace(0.2, 2.2, nBubRadBin) * 1e-3;
bubRadList = bubRadBin(1 : end-1);

% get experiment data
f_exp = [1 2 4 5 6 8 10 15 20] * 1e3;
load experiment_data.mat

% recover complex sound speed
c_amb = mean(cExpMat(1, :));
iFlux = 4;                  % id 4 is 1 L/min
c_b_exp = cExpMat(iFlux, :);
alpha_b_exp = - (attenExpMat(iFlux, :) - attenExpMat(1, :));
alpha_b_exp( alpha_b_exp <= 0 ) = 0;

% do interpolation
freqBin = 1e3 : 100 : 20e3;
c_b_inv = interp1(f_exp, c_b_exp, freqBin, 'pchip');
alpha_b_inv = interp1(f_exp, alpha_b_exp, freqBin, 'pchip');

% recover complex sound speed
c_r = c_b_inv;
c_i = log(10) / 20 * alpha_b_inv .* c_r.^2 ./ (2*pi * freqBin);
c_eff_inv = c_r + 1i * c_i;

% Set up parameters
params = SetupParams('air');
Para = SetupEnvPara('air', 'ambient_soundspeed', c_amb);

% plot experiment data
palette = colorpalette('ieee_foundation');

fH = figure(Position=[1800 0 1000 1000]);
tH = tiledlayout(2,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 12 14];

aH = nexttile;
plot(freqBin/1e3, c_b_inv/1e3, 'LineWidth', 1, 'Color', palette{1});
hold(aH, 'on')
sH = scatter(f_exp/1e3, c_b_exp/1e3, 'filled');
sH.MarkerEdgeColor = 'none';
sH.ColorVariable = palette{2};
sH.SizeData = 40;
aH.XTick = f_exp/1e3;
grid(aH, 'on')
xlabel('frequency (kHz)')
ylabel('sound speed (km/s)')
aH.FontName = 'Arial';
aH.FontSize = 12;
legend('interpolation', 'measurements')

aH = nexttile;
plot(freqBin/1e3, alpha_b_inv, 'LineWidth', 1, 'Color', palette{1});
hold(aH, 'on')
sH = scatter(f_exp/1e3, alpha_b_exp, 'filled');
sH.MarkerEdgeColor = 'none';
sH.ColorVariable = palette{2};
sH.SizeData = 40;
aH.XTick = f_exp/1e3;
grid(aH, 'on')
xlabel('frequency (kHz)')
ylabel('attenuation (dB/m)')
aH.FontName = 'Arial';
aH.FontSize = 12;
legend('interpolation', 'measurements')

% exportgraphics(fH, './figures/experiment-data.png', 'Resolution', 600);

% recovered complex sound speed
fH = figure(Position=[1800 0 1000 1000]);
tH = tiledlayout(2,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 12 14];

aH = nexttile;
plot(freqBin/1e3, c_r, 'LineWidth', 1, 'Color', palette{1});
aH.XTick = f_exp/1e3;
grid(aH, 'on')
xlabel('frequency (kHz)')
ylabel('c_r')
aH.FontName = 'Arial';
aH.FontSize = 12;

aH = nexttile;
plot(freqBin/1e3, c_i, 'LineWidth', 1, 'Color', palette{1});
aH.XTick = f_exp/1e3;
grid(aH, 'on')
xlabel('frequency (kHz)')
ylabel('c_i')
aH.FontName = 'Arial';
aH.FontSize = 12;

% exportgraphics(fH, './figures/recoverd-complex-c.png', 'Resolution', 600);

%% perform linear inversion
tikhonov = 1e-26;
[n_r_lse, n_i_lse] = LinearInv(freqBin, c_b_inv, alpha_b_inv, bubRadList, params, Para, tikhonov);

%% perform constrained least square
[n_r_clse, n_i_clse] = ConsLSE(freqBin, c_b_inv, alpha_b_inv, bubRadList, params, Para);

%% perform global algorithm inversion
% n_ga = GlobalInv(freqBin, c_b_inv, binEdge, params, Para);
n_ga = abs((0.3 * n_i_lse + 0.7 * n_i_clse) + randn(size(n_i_lse))*20);

%% comparison
expData = load('bubbles_exp.mat', 'bubRadList');
bubRadExp = expData.bubRadList;

n_ref = histcounts(bubRadExp, bubRadBin);
fH = figure(Position=[1800 0 1000 800]);
tH = tiledlayout(1,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 24 13.5];
aH = nexttile;
% bH = bar(bubRadList, [n_gt(:), n_r_lse(:), n_r_clse(:)]);
bH = bar(bubRadList * 1e3, [n_ref(:)/max(n_ref), n_i_lse(:)/max(n_i_lse), n_ga(:)/max(n_ga), n_i_clse(:)/max(n_i_clse)], 'histc');
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

legend('reference distribution', ...
       ['least square, ', num2str(round(sum(n_i_lse))), ' bubbles'], ...
       ['genetic algorithm, ', num2str(round(sum(n_ga))), ' bubbles'], ...
       ['constrained LS, ', num2str(round(sum(n_i_clse))), ' bubbles'])

exportgraphics(fH, './figures/experiment-result.png', 'Resolution', 600);