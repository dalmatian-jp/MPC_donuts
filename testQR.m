% Define selected Q and R values
Q = diag([6457.5, 7801.7, 7838.3, 8860.4]);  % Now Q is 4×4
R = diag([0.02, 0.01]);  % Now R is 2×2

% Initial state and reference state
x = [0.0873; 0; 0; 0];  % Current state
x_ref = [0; 0; 0; 0];  % Desired equilibrium state

% Time step and prediction horizon
N = 10;  % Prediction horizon
dt = 0.01; % Sampling time (s)

% Call the MPC controller with selected Q and R
u_opt = mpcController(x, x_ref, systemParams, N, Q, R, dt);

% Display the optimal control input
disp('Optimal control input:');
disp(u_opt);

% Simulate the next state
x_next = doubleLinkDynamics(0, x, u_opt, systemParams);

% Display the next state
disp('Next state:');
disp(x_next);