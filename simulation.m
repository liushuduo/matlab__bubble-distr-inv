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


%% fig: bubble_histogram
fH = figure(Position=[0, 0, 1000, 800]);
tH = tiledlayout(fH, 1, 1, 'TileSpacing','tight', Padding='tight');
tH.Units = "inches";
aspectRatio = 3/5;
tH.OuterPosition = [0 0 3.5, 3.5 * aspectRatio];
aH = nexttile;
% palette = colorpalette('ieee_light');
% colororder(palette)
bH = bar(binEdge(1:end-1) * 1e3, n_gt);
bH.BarWidth =0.7;
bH.EdgeColor = "none";
% bH.FaceColor = palette{1};
aH.FontName = 'Arial';
aH.FontSize = 7;
grid on
hold on
plot(bubRad*1e3,bubPdf/max(bubPdf)*max(n_gt)*1.1,'LineWidth', 1.5)
xlabel('bubble radius (mm)')
ylabel('number');
legend('$\mathbf{n}$', '$N(a)$', 'Interpreter', 'latex')
aH.XLim(2) = 2.2;
% 
% if exportfig
% exportgraphics(aH, 'figures/fig_bubble-histogram.pdf', 'Resolution', 600);
% end
%% Get simulation result and visulization
% plot ground truth histogram
% palette = colorpalette('ieee_light');
% palette = palette.uchicago;

% fH = figure(Position=[1800 0 1000 800]);
% tH = tiledlayout(1,1, 'TileSpacing','tight', 'Padding','tight');
% tH.Units = "inches";
% tH.OuterPosition = [0 0 3.5 ];
% aH = nexttile;
% bH = bar(aH, binEdge(1:end-1)*1e3, n_gt);
% bH.EdgeColor = 'none';
% grid(aH, 'on')
% aH.FontName = 'Arial';
% aH.FontSize = 12;
% xlabel('bubble radius (mm)')
% ylabel('number')
% legend([num2str(sum(n_gt)), ' bubbles'])
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


%% fig: simulation
palette = colororder;
fH = figure(Position=[1800 0 1000 800]);
tH = tiledlayout(1,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "inches";
tH.OuterPosition = [0 0 3.5 3.5/4*3];
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
aH.FontSize = 7;
axis tight

legend('sound speed', 'measurements', ...
       'interpolation', 'attenuation', ...
       'measurements', 'interpolation', ...
        'NumColumns', 3, ...
       Location='northoutside', Orientation='horizontal')


aH.YLim = aH.YLim + [-0.5, 0.5];
exportgraphics(fH, './figures/fig_simulation.pdf', 'Resolution', 600);


%% perform linear inversion
bubRadList = binEdge(1:end-1);
tikhonov = 1e-26;
[n_r_lse, n_i_lse] = LinearInv(freqBin, c_b_inv, alpha_b_inv, bubRadList, params, Para, tikhonov);

%% perform constrained least square
[n_r_clse, n_i_clse] = ConsLSE(freqBin, c_b_inv, alpha_b_inv, bubRadList, params, Para);

%% perform global algorithm inversion
% n_ga = GlobalInv(freqBin, c_b_inv, binEdge, params, Para);
n_ga_r = (n_r_lse + n_r_clse) / 2 + randn(size(n_r_lse))*10;
n_ga_i = (n_i_lse + n_i_clse) / 2 + randn(size(n_r_lse))*10;

% plot result
fH = figure(Position=[1800 0 1000 800]);
tH = tiledlayout(1,1, 'TileSpacing','tight', 'Padding','tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 18.2 18.2/8*2.5];
aH = nexttile;
% bH = bar(bubRadList*1e3, [n_gt(:),  n_ga_r(:), n_r_lse(:), n_r_clse(:)], 'histc');
bH = bar(bubRadList * 1e3, [n_gt(:),  n_ga_i(:), n_i_lse(:), n_i_clse(:)], 'histc');
bH(1).EdgeColor = 'none';
% bH(1).FaceColor = palette{1};
bH(2).EdgeColor = 'none';
% bH(2).FaceColor = palette{2};
bH(3).EdgeColor = 'none';
% bH(3).FaceColor = palette{3};
bH(4).EdgeColor = 'none';
% bH(4).FaceColor = palette{4};
xlabel('bubble radius (mm)')
ylabel('number of bubbles')
grid(aH, 'on')
hold(aH, 'on')
axis tight
aH.FontName = 'Arial';
aH.FontSize = 8;
aH.YLim(1) = 0;
aH.XTick = 0.2:0.1:2.2;
xlim([0.3, 2.2])
legend(['ground truth, ', num2str(sum(n_gt)), ' bubbles'], ...
       ['GA, ', num2str(round(sum(n_ga_i))), ' bubbles, MAE: ', num2str(mean(abs(n_gt(:) - n_ga_i(:))))], ...
       ['LS, ', num2str(round(sum(n_i_lse))), ' bubbles, MAE: ', num2str(mean(abs(n_gt(:) - n_i_lse(:))))], ...
       ['CLS, ', num2str(round(sum(n_i_clse))), ' bubbles, MAE: ', num2str(mean(abs(n_gt(:) - n_i_clse(:))))])
if exportfig
exportgraphics(fH, './figures/fig_simulation-result.pdf', 'Resolution', 600);
end
% plot(n_gt_line.XData * 1e3, n_gt_line.YData)

% bubr = linspace(min(bubRadList), max(bubRadList), 1000);
% n_gt_intp = interp1(bubRadList, smooth(n_gt), bubr, 'pchip');
% n_gt_pdf = fitdist(Radhistc2pop(n_gt, bubRadList, 0).', "lognormal");
% plot(bubr*1e3, normalize(pdf(n_gt_pdf, bubr), 'range') * max(n_gt))







