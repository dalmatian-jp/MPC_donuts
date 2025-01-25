systemParams = struct(...
    'm1',   25.00, ...
    'm2',   50.00, ...
    'L1',   0.90, ...
    'L2',   0.85, ...
    'I1',   1.00, ...
    'I2',   1.20, ...
    'com1', 0.45, ...
    'com2', 0.42, ...
    'g',    9.81);

% データの抽出
data = readtable('EKF_Simulation_AllData.csv');

q1_s = data.q1_s_true;
q2_s = data.q2_s_true;
dq1_s = data.dq1_s_true;
dq2_s = data.dq2_s_true;
u1_s = data.u1_s;
u2_s = data.u2_s;

q1_l = data.q1_l_true;
q2_l = data.q2_l_true;
dq1_l = data.dq1_l_true;
dq2_l = data.dq2_l_true;
u1_l = data.u1_l;
u2_l = data.u2_l;

% パラメータの設定（pvstate）
pvstate = [0.90; 0.45; 25.00; 1.00; 0.85; 0.42; 50.00; 1.20; 9.81];

% COPを格納する配列
num_steps = height(data);
COP_s_values = zeros(num_steps, 1);
COP_l_values = zeros(num_steps, 1);

% 各ステップでCOPを計算
for i = 1:num_steps
    % 状態ベクトルとトルクベクトルの作成
    x_s = [q1_s(i); q2_s(i); dq1_s(i); dq2_s(i)];
    tau_s = [u1_s(i); u2_s(i)];

    x_l = [q1_l(i); q2_l(i); dq1_l(i); dq2_l(i)];
    tau_l = [u1_l(i); u2_l(i)];
    
    % COPの計算（makeCOP関数を使用）
    COP_s_values(i) = makeCOP(x_s, tau_s, systemParams, pvstate);
    COP_l_values(i) = makeCOP(x_l, tau_l, systemParams, pvstate);
end

% 結果の表示
figure;
plot(COP_s_values, 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('COP');
title('Center of Pressure (COP) Over Time');
grid on;

data_wt_cop = addvars(data, COP_s_values, COP_l_values, ...
               'NewVariableNames', {'COP_s', 'COP_l'});

% CSVファイルにヘッダー付きで保存
writetable(data_wt_cop, 'EKF_Simulation_AllData_with_COP.csv');
