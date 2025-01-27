clc; close all;
[estSmall, s_StateHistory, s_uHistory, estLarge, l_StateHistory, l_uHistory] = twoPerturbationsEKF(6000);
allDataMatrix = [estSmall, s_StateHistory, s_uHistory, ...
                 estLarge, l_StateHistory, l_uHistory];

% CSVファイルに書き出し
writematrix(allDataMatrix, 'data/AllData_takami_60s.csv');

disp('すべてのシミュレーションデータを1つのCSVファイルに保存しました。');
