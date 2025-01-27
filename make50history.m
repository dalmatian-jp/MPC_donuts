clc; close all;

for i = 1:50
    [estSmall, s_StateHistory, s_uHistory] = for50(3000);
    allDataMatrix = [estSmall, s_StateHistory, s_uHistory];
                    %  estLarge, l_StateHistory, l_uHistory];

    % CSVファイルに書き出し
    filename = sprintf('new50\\EKF_Simulation_AllData_%02d.csv', i);
    writematrix(allDataMatrix, filename);

    disp(['シミュレーションデータを ' filename ' に保存しました。']);
end
