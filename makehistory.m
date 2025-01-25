clc; close all;
[estSmall, s_StateHistory, s_uHistory, estLarge, l_StateHistory, l_uHistory] = twoPerturbationsEKF(3000);
allDataMatrix = [estSmall, s_StateHistory, s_uHistory, ...
                 estLarge, l_StateHistory, l_uHistory];

% CSVファイルに書き出し
writematrix(allDataMatrix, 'EKF_Simulation_AllData.csv');

disp('すべてのシミュレーションデータを1つのCSVファイルに保存しました。');
