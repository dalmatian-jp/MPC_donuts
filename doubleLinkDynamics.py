# Re-import necessary modules after environment reset
import pandas as pd
import numpy as np
import os

# Ensure the directory exists
os.makedirs('induced_anal', exist_ok=True)

# Load the simulation data
file_path = 'data/EKF_Simulation_AllData.csv'
data = pd.read_csv(file_path)

# Define the double-link dynamics function again
def double_link_dynamics(t, x, u, params):
    m1 = params['m1']
    m2 = params['m2']
    L1 = params['L1']
    L2 = params['L2']
    I1 = params['I1']
    I2 = params['I2']
    g = params['g']

    q1, q2, dq1, dq2 = x

    # Mass matrix (M)
    M11 = m1 * L1**2 / 3 + m2 * (L1**2 + (L2**2) / 3 + L1 * L2 * np.cos(q2)) + I1 + I2
    M12 = m2 * ((L2**2) / 3 + L1 * L2 * np.cos(q2)) + I2
    M21 = M12
    M22 = m2 * (L2**2) / 3 + I2
    M = np.array([[M11, M12], [M21, M22]])

    # Coriolis and centrifugal terms (C)
    C1 = -m2 * L1 * L2 * np.sin(q2) * dq2**2 - 2 * m2 * L1 * (L2 / 2) * np.sin(q2) * dq1 * dq2
    C2 = m2 * L1 * L2 * np.sin(q2) * dq1**2
    C = np.array([C1, C2])

    # Gravitational terms (F)
    F1 = (m1 * L1 / 2 + m2 * L1) * g * np.sin(q1) + m2 * (L2 / 2) * g * np.sin(q1 + q2)
    F2 = m2 * (L2 / 2) * g * np.sin(q1 + q2)
    F = np.array([F1, F2])

    # Compute angular accelerations (ddq)
    ddq = np.linalg.inv(M + 1e-6 * np.eye(2)) @ (u - C - F)

    return np.array([dq1, dq2, ddq[0], ddq[1]])

# Example parameters for the double inverted pendulum
params_example = {
    'm1': 1.0,     # Mass of link 1 (kg)
    'm2': 1.0,     # Mass of link 2 (kg)
    'L1': 1.0,     # Length of link 1 (m)
    'L2': 1.0,     # Length of link 2 (m)
    'I1': 0.1,     # Moment of inertia of link 1 (kg·m²)
    'I2': 0.1,     # Moment of inertia of link 2 (kg·m²)
    'g': 9.81      # Gravitational acceleration (m/s²)
}

# Extract relevant columns for small and large perturbations
small_perturbation = data[["q1_s_true", "q2_s_true", "dq1_s_true", "dq2_s_true", "u1_s", "u2_s"]]
large_perturbation = data[["q1_l_true", "q2_l_true", "dq1_l_true", "dq2_l_true", "u1_l", "u2_l"]]

# Analyze the small perturbation data to compute induced accelerations
induced_accelerations_small = []
for _, row in small_perturbation.iterrows():
    x = np.array([row["q1_s_true"], row["q2_s_true"], row["dq1_s_true"], row["dq2_s_true"]])
    u = np.array([row["u1_s"], row["u2_s"]])
    dxdt = double_link_dynamics(0, x, u, params_example)
    induced_accelerations_small.append(dxdt)

induced_accelerations_small_df = pd.DataFrame(induced_accelerations_small, columns=["dq1", "dq2", "ddq1", "ddq2"])

# Analyze the large perturbation data similarly
induced_accelerations_large = []
for _, row in large_perturbation.iterrows():
    x = np.array([row["q1_l_true"], row["q2_l_true"], row["dq1_l_true"], row["dq2_l_true"]])
    u = np.array([row["u1_l"], row["u2_l"]])
    dxdt = double_link_dynamics(0, x, u, params_example)
    induced_accelerations_large.append(dxdt)

induced_accelerations_large_df = pd.DataFrame(induced_accelerations_large, columns=["dq1", "dq2", "ddq1", "ddq2"])

# Display the results to the user
print("Induced Accelerations Analysis (Small Perturbations)")
print(induced_accelerations_small_df)
print("\nInduced Accelerations Analysis (Large Perturbations)")
print(induced_accelerations_large_df)

induced_accelerations_small_df.to_csv('induced_anal/Induced_Accelerations_Small.csv', index=False)
induced_accelerations_large_df.to_csv('induced_anal/Induced_Accelerations_Large.csv', index=False)
print("Data saved to 'induced_anal/Induced_Accelerations_Small.csv' and 'induced_anal/Induced_Accelerations_Large.csv'")


