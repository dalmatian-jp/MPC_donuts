import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def makeCOP(x, u, params, pvstate):
    q1, q2, dq1, dq2 = x

    l1, s1, m1, J1, l2, s2, m2, J2, g = pvstate
    q0 = 4 * np.pi / 180  # Convert degrees to radians
    l0 = 0.27
    k1 = 0.468
    k2 = 0.73

    dxdt = double_link_dynamics_takami(0, x, u, params)  # Use existing dynamics function
    ddq1 = dxdt[2]  # Extract angular acceleration of q1
    ddq2 = dxdt[3]  # Extract angular acceleration of q2

    ut = ((-(m1 * k1**2 + m2 * l1**2 + m2 * k2**2) + 2 * m2 * l1 * k2 * np.cos(q2) +
          (m1 * l0 * k1 + m2 * l0 * l1) * np.cos(q0) * np.cos(q1 - q0) -
          m2 * l0 * k2 * np.cos(q0) * np.cos(q2 + q1 - q0)) * ddq1 +
         (-m2 * k2**2 + m2 * l1 * k2 * np.cos(q2) - 
          m2 * l0 * k2 * np.cos(q0) * np.cos(q2 + q1 - q0)) * ddq2)

    fv = ((m1 * k1 + m2 * l1) * np.cos(q1 - q0) * ddq1 -
         (m2 * k2 + m2 * l2) * np.cos(q2 + q1 - q0) * (ddq2 + ddq1) + (m1 + m2) * g)

    cop = ut / fv if fv != 0 else np.nan  # Avoid division by zero

    return cop


def double_link_dynamics_takami(t, x, u, params):
    m1, m2 = params['m1'], params['m2']
    L1, L2 = params['L1'], params['L2']
    I1, I2 = params['I1'], params['I2']
    g = params['g']

    r1 = 0.64 * L1
    r2 = 0.36 * L2
    q1, q2, dq1, dq2 = x
    q1_new = q1 + np.pi / 2

    a = I1 + I2 + m1 * r1**2 + m2 * (L1**2 + r2**2)
    b = m2 * L1 * r2
    d = I2 + m2 * r2**2

    s1 = np.sin(q1_new)
    s2 = np.sin(q2)
    c1 = np.cos(q1_new)
    c2 = np.cos(q2)
    c12 = np.cos(q1_new + q2)

    M11 = a + 2 * b * c2
    M12 = d + b * c2
    M = np.array([[M11, M12],
                  [M12, d]])

    C = np.array([
        [-b * s2 * dq2, -b * s2 * (dq1 + dq2)],
        [b * s2 * dq2, 0]
    ])

    G1 = -g * ((m1 * r1 + m2 * L1) * c1 + m2 * r2 * c12)
    G2 = -g * m2 * r2 * c12
    G = np.array([G1, G2])

    ddq = np.linalg.solve(M + 1e-6 * np.eye(2), u - C @ np.array([dq1, dq2]) - G)
    dxdt = np.array([dq1, dq2, ddq[0], ddq[1]])

    return dxdt

systemParams = {
    'm1': 28.00, 'm2': 53.00,
    'L1': 0.90, 'L2': 0.88,
    'I1': 9.21, 'I2': 5.35,
    'com1': 0.58, 'com2': 0.32,
    'g': 9.81
}
pvstate = [0.90, 0.58, 28.00, 9.21, 0.88, 0.32, 53.00, 5.35, 9.81]


# システムパラメータの定義
systemParams = {
    'm1': 28.00, 'm2': 53.00,
    'L1': 0.90, 'L2': 0.88,
    'I1': 9.21, 'I2': 5.35,
    'com1': 0.58, 'com2': 0.32,
    'g': 9.81
}

# pvstate の設定
pvstate = [0.90, 0.58, 28.00, 9.21, 0.88, 0.32, 53.00, 5.35, 9.81]

# データが格納されているフォルダーのパス
folder_path = "new50"
output_folder = "processed_data"  # 結果を保存するフォルダー
os.makedirs(output_folder, exist_ok=True)  # フォルダーがなければ作成

# CSVファイルの一覧を取得
csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

# 処理する列のリスト
columns_to_convert = [
    'q1_s_true', 'q2_s_true', 'dq1_s_true', 'dq2_s_true',
    'u1_s', 'u2_s'
]

# 各CSVファイルに対してCOPを計算
for file_name in csv_files:
    file_path = os.path.join(folder_path, file_name)

    # データの読み込み
    data = pd.read_csv(file_path)

    # 必要な列を数値型に変換
    for col in columns_to_convert:
        if col in data.columns:
            data[col] = pd.to_numeric(data[col], errors='coerce')

    # NaNのチェック
    if data.isnull().values.any():
        print(f"Warning: NaN values detected in {file_name}. Check your data!")

    # 各変数の抽出
    q1_s = data['q1_s_true'].values
    q2_s = data['q2_s_true'].values
    dq1_s = data['dq1_s_true'].values
    dq2_s = data['dq2_s_true'].values
    u1_s = data['u1_s'].values
    u2_s = data['u2_s'].values

    # COPを格納する配列
    num_steps = len(data)
    COP_s_values = np.zeros(num_steps)

    # 各ステップでCOPを計算
    for i in range(num_steps):
        x_s = [q1_s[i], q2_s[i], dq1_s[i], dq2_s[i]]
        tau_s = [u1_s[i], u2_s[i]]

        # COPの計算（makeCOP関数を使用）
        COP_s_values[i] = makeCOP(x_s, tau_s, systemParams, pvstate)

    # データにCOPの列を追加
    data['COP_s'] = COP_s_values

    # 結果の保存パスを設定
    output_file_path = os.path.join(output_folder, file_name)

    # CSVファイルに保存
    data.to_csv(output_file_path, index=False)

    print(f"Processed: {file_name} → Saved to {output_file_path}")

# 結果のプロット（最後のCSVファイルを例として）
plt.figure(figsize=(10, 5))
plt.plot(COP_s_values, label="COP_s", linewidth=1.5)
plt.xlabel("Time Step")
plt.ylabel("COP")
plt.title(f"Center of Pressure (COP) for {file_name}")
plt.legend()
plt.grid(True)
plt.show()
