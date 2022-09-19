%% from N(a) to N
clear; close all; clc
palette = colorpalette("ieee_foundation");

% the N(a) is a function where only one or zero exists
load bubbles_exp.mat
fH = figure(Position=[1800 -200 1000 500]);
tH = tiledlayout(1,1, 'TileSpacing','tight', Padding='tight');
tH.Units = "inches";
tH.OuterPosition = [0 0 3.5 3.5/16 * 9];
aH = nexttile;
sH = stem(aH, bubRadList*1e3, ones(size(bubRadList)));
sH.Color = palette{1};
sH.Marker = "none";
xlabel('bubble radius (mm)')
ylabel('$N(a)$', 'Interpreter','latex')
ylim([0 1.5])
grid on
aH.FontName = 'Times';
aH.FontSize = 8;
exportgraphics(fH, './figures/Na.png', 'Resolution', 600);

% convert N(a) to vector N
nBin = 31;
binEdge = linspace(min(bubRadList), max(bubRadList), nBin);
N = histcounts(bubRadList, binEdge);

% plot the histogram
fH = figure(Position=[1800 -200 1000 500]);
tH = tiledlayout(1,1, 'TileSpacing','tight', Padding='tight');
tH.Units = "inches";
tH.OuterPosition = [0 0 3.5 3.5/16*9];
aH = nexttile;
bH = bar(binEdge(1:end-1) * 1e3, N);
aH.YAxis.Exponent = 3;
bH.EdgeColor = "none";
bH.FaceColor = palette{1};
aH.FontName = 'Times';
aH.FontSize = 8;
grid on
xlabel('bubble radius (mm)')
ylabel('$\mathbf{n}$', 'Interpreter','latex');
exportgraphics(fH, './figures/vecN.png', 'Resolution', 600);



