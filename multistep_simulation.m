% Number of steps to simulate
numSteps = 400;
dt       = 0.01;   % Example sampling time; ensure it matches your MPC
x        = [0.0873; 0; 0; 0];  % Initial state
x_ref    = [0; 0; 0; 0];       % Reference state

Q = diag([6457.5, 7801.7, 7838.3, 8860.4]);  % Q is 4×4
R = diag([0.02, 0.01]);                      % R is 2×2

% Initialize storage for state history
stateHistory = zeros(numSteps + 1, length(x));
stateHistory(1, :) = x';

% Simulate over multiple steps
for k = 1:numSteps
    
    % 1) Compute optimal control input
    u_opt = mpcController(x, x_ref, systemParams, N, Q, R, dt);

    % 2) Compute the state derivative at the current instant
    dxdt = doubleLinkDynamics(0, x, u_opt, systemParams);

    % 3) Update the state via simple Euler integration
    x = x + dt * dxdt;

    % Store the updated state
    stateHistory(k + 1, :) = x';
end

% Plot state evolution
figure;
plot(0:numSteps, stateHistory(:, 1), '-o', 'DisplayName', 'q1');
hold on;
plot(0:numSteps, stateHistory(:, 2), '-x', 'DisplayName', 'q2');
plot(0:numSteps, stateHistory(:, 3), '-s', 'DisplayName', 'q1dot');
plot(0:numSteps, stateHistory(:, 4), '-d', 'DisplayName', 'q2dot');

% Also plot the reference lines for q1 & q2 (and optionally q1dot, q2dot if reference is nonzero)
plot(0:numSteps, repmat(x_ref(1), numSteps + 1, 1), '--', 'DisplayName', 'q1 Ref');
plot(0:numSteps, repmat(x_ref(2), numSteps + 1, 1), '--', 'DisplayName', 'q2 Ref');

xlabel('Time Step');
ylabel('State Value');
legend();
title('State Convergence to Reference');
grid on;
