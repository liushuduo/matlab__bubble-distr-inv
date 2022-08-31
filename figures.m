%% from N(a) to N
clear; close all; clc
palette = GetPalette();
palette = palette.uchicago;

% the N(a) is a function where only one or zero exists
load bubbles_exp.mat
fH = figure(Position=[1800 -200 1000 500]);
tH = tiledlayout(1,1, 'TileSpacing','tight', Padding='tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 8.6 5.4] * 0.95;
aH = nexttile;
sH = stem(aH, bubRadList*1e3, ones(size(bubRadList)));
sH.Color = palette{5};
sH.Marker = "none";
xlabel('bubble radius')
ylabel('$N(a)$', 'Interpreter','latex')
ylim([0 1.5])
aH.FontName = 'Arial';
aH.FontSize = 12;
aH.YTick = [];   aH.XTick = [];
% exportgraphics(fH, './figures/Na.png', 'Resolution', 600);

% convert N(a) to vector N
nBin = 31;
binEdge = linspace(min(bubRadList), max(bubRadList), nBin);
N = histcounts(bubRadList, binEdge);

% plot the histogram
fH = figure(Position=[1800 -200 1000 500]);
tH = tiledlayout(1,1, 'TileSpacing','tight', Padding='tight');
tH.Units = "centimeters";
tH.OuterPosition = [0 0 8.6 5.4] * 0.95;
aH = nexttile;
bH = bar(binEdge(1:end-1), N);
bH.EdgeColor = "none";
bH.FaceColor = palette{5};
aH.FontName = 'Arial';
aH.FontSize = 12;
aH.YTick = [];   aH.XTick = [];
xlabel('bubble radius')
ylabel('$\mathbf{n}$', 'Interpreter','latex');
exportgraphics(fH, './figures/vecN.png', 'Resolution', 600);



