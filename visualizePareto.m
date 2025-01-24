% Assumes 'optimalParams' and 'optimalObjectives' are already available
% OptimalParams: [35 x 6] matrix containing Q and R values
% OptimalObjectives: [35 x 2] matrix containing RMSE and Energy values

% Extract Q and R values
Q_values = optimalParams(:, 1:4); % First 4 columns are Q
R_values = optimalParams(:, 5:6); % Last 2 columns are R
RMSE = optimalObjectives(:, 1);   % First column of objectives
Energy = optimalObjectives(:, 2); % Second column of objectives

% Pareto Front Scatter Plot
figure;
scatter(RMSE, Energy, 100, 'filled');
xlabel('RMSE');
ylabel('Energy');
title('Pareto Front: RMSE vs. Energy');
grid on;
colorbar;
caxis([min(Energy), max(Energy)]);
text(RMSE, Energy, string(1:length(RMSE)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Heatmap of Q and R Values
figure;
subplot(1, 2, 1);
imagesc(Q_values);
title('Q Values (Heatmap)');
xlabel('Q Components');
ylabel('Solutions');
colorbar;
xticks(1:4);
xticklabels({'Q1', 'Q2', 'Q3', 'Q4'});

subplot(1, 2, 2);
imagesc(R_values);
title('R Values (Heatmap)');
xlabel('R Components');
ylabel('Solutions');
colorbar;
xticks(1:2);
xticklabels({'R1', 'R2'});

% RMSE and Energy Histograms
figure;
subplot(1, 2, 1);
histogram(RMSE, 10, 'FaceColor', 'b');
xlabel('RMSE');
ylabel('Frequency');
title('RMSE Distribution');
grid on;

subplot(1, 2, 2);
histogram(Energy, 10, 'FaceColor', 'r');
xlabel('Energy');
ylabel('Frequency');
title('Energy Distribution');
grid on;
