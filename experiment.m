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

%% fig:exp
fH = figure(Position=[1800 0 1000 800]);
tH = tiledlayout(1,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "inches";
tH.OuterPosition = [0 0 3.5 3.5/4*3];
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
aH.FontSize = 7;
axis tight

% legend('sound speed measurements','$\mathbf{c}_b$', ...
%     'attenuation measurements', '$\boldsymbol{\alpha}_b$', 'Interpreter', 'latex', ...
%     'NumColumns', 2, ...
%        Location='northoutside', Orientation='horizontal')

legend('sound speed measurements', 'interpolation', ...
       'attenuation measurements', 'interpolation', ...
        'NumColumns', 2, ...
       Location='northoutside', Orientation='horizontal')

xlim([0, 9])
aH.YLim = aH.YLim + [-0.5, 0.5];
exportgraphics(fH, './figures/fig_exp.pdf', 'Resolution', 600);


%% Set up parameters
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
% fH = figure(Position=[1800 0 1000 1000]);
% tH = tiledlayout(2,1, 'TileSpacing','tight', 'Padding','tight');
% tH.Units = "centimeters";
% tH.OuterPosition = [0 0 12 14];
% 
% aH = nexttile;
% plot(freqBin/1e3, c_r, 'LineWidth', 1, 'Color', palette{1});
% aH.XTick = f_exp/1e3;
% grid(aH, 'on')
% xlabel('frequency (kHz)')
% ylabel('c_r')
% aH.FontName = 'Arial';
% aH.FontSize = 12;
% 
% aH = nexttile;
% plot(freqBin/1e3, c_i, 'LineWidth', 1, 'Color', palette{1});
% aH.XTick = f_exp/1e3;
% grid(aH, 'on')
% xlabel('frequency (kHz)')
% ylabel('c_i')
% aH.FontName = 'Arial';
% aH.FontSize = 12;

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

% All bubble radii for current simulation
load('n_gt.mat')
n_ref = n_gt;

fH = figure(Position=[1800 0 1000 800]);
tH = tiledlayout(1,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 18.2 18.2/8*2.5];
aH = nexttile;
% bH = bar(bubRadList, [n_gt(:), n_r_lse(:), n_r_clse(:)]);
bH = bar(bubRadList * 1e3, [n_ref(:)/max(n_ref),n_ga(:)/max(n_ga),  n_i_lse(:)/max(n_i_lse), n_i_clse(:)/max(n_i_clse)], 'histc');
% bH = bar(bubRadList, [n_gt(:), n_r_lse(:), n_r_clse(:)]);
bH(1).EdgeColor = 'none';
% bH(1).FaceColor = palette{1};
bH(2).EdgeColor = 'none';
% bH(2).FaceColor = palette{2};
bH(3).EdgeColor = 'none';
% bH(3).FaceColor = palette{3};
bH(4).EdgeColor = 'none';
% bH(4).FaceColor = palette{4};
xlabel('bubble radius (mm)')
ylabel('normalized number of bubbles')
grid(aH, 'on')
hold(aH, 'on')
axis tight
aH.FontName = 'Arial';
aH.FontSize = 8;
aH.YLim(1) = 0;
xlim([0.3, 2.2])
aH.XTick = 0.4:0.1:2.2;

legend('reference distribution', ...
       ['GA, ', num2str(round(sum(n_ga))), ' bubbles'], ...
       ['LS, ', num2str(round(sum(n_i_lse))), ' bubbles'], ...
       ['CLS, ', num2str(round(sum(n_i_clse))), ' bubbles'])

exportgraphics(fH, './figures/fig_exp-result.pdf', 'Resolution', 600);