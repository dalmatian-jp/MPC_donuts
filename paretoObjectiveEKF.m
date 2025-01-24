function objectives = paretoObjectiveEKF(params, x0, x_ref, N, dt, systemParams)
% PARETOOBJECTIVEEKF:
%   Objective function for multi-objective GA that runs an EKF + MPC closed-loop
%   for a given Q,R. Outputs [RMSE, Energy].
%
% Inputs:
%   params  = [q1, q2, q3, q4, r1, r2] (6x1 vector)
%   x0      = initial state (4x1)
%   x_ref   = reference state (4x1)
%   N       = MPC horizon
%   dt      = sampling time
%   systemParams = struct with dynamics parameters
%
% Output:
%   objectives = [ RMSE, totalEnergy ]

    % 1) Extract Q,R from 'params'
    Q = diag(params(1:4));  % e.g. Q for the 4 states
    R = diag(params(5:6));  % R for the 2 inputs

    % 2) Define how many simulation steps
    numSteps = 200;  % or however long you want to test each candidate

    % 3) Define/Initialize EKF stuff
    Q_ekf = 1e-6 * eye(4);   % or your chosen process noise
    R_ekf = diag([1e-5, 1e-5, 1e-4, 1e-4]); % measurement noise
    x_true = x0;
    x_hat  = x_true + 0.01*randn(4,1);
    P      = 1e-3*eye(4);   % initial covariance of EKF

    % 4) Tracking variables
    totalError  = 0;
    totalEnergy = 0;

    for k = 1:numSteps
        % (A) Generate a noisy measurement from the true state
        z_meas = x_true + mvnrnd(zeros(1,4), R_ekf)';

        % (B) EKF Predict step
        u_pred = mpcController(x_hat, x_ref, systemParams, N, Q, R, dt);
        [xPred, PPred] = ekfPredict(x_hat, P, ...
            @(xx,uu) fEuler(xx,uu,systemParams,dt), ...
            @(xx,uu) A_jacobian(xx,uu,systemParams,dt), ...
            zeros(4,1), Q_ekf, u_pred);

        % (C) EKF Update step
        [x_hat, P] = ekfUpdate(xPred, PPred, z_meas, ...
            @(xx) h_meas(xx), @(xx) H_jacobian(xx), R_ekf);

        % (D) Use the updated estimate in the real control
        u_opt = mpcController(x_hat, x_ref, systemParams, N, Q, R, dt);

        % (E) Propagate TRUE system
        dxdt_true = doubleLinkDynamics(0, x_true, u_opt, systemParams);
        x_true = x_true + dt*dxdt_true + mvnrnd(zeros(1,4), Q_ekf)';

        % (F) Accumulate error (true state vs. reference)
        totalError  = totalError + norm(x_true - x_ref)^2;
        % (G) Accumulate energy usage
        totalEnergy = totalEnergy + norm(u_opt)^2;
    end

    % Final objectives
    RMSE   = sqrt( totalError / numSteps );
    Energy = totalEnergy;

    objectives = [RMSE, Energy];
end


%% ========================================================================
% The discrete update used in the EKF:
function x_next = fEuler(x, u, params, dt)
    dxdt = doubleLinkDynamics(0, x, u, params);
    x_next = x + dt*dxdt;
end

function A = A_jacobian(x, u, params, dt)
    syms q1 q2 dq1 dq2 tau1 tau2 real
    state = [q1; q2; dq1; dq2];
    ctrl  = [tau1; tau2];

    dxdt_sym = doubleLinkDynamics(0, state, ctrl, params);
    A_sym    = jacobian(dxdt_sym, state);

    A_cont = double(subs(A_sym, [state; ctrl], [x; u]));
    A      = eye(4) + dt*A_cont;
end

function z = h_meas(x)
    z = x;
end

function H = H_jacobian(~)
    H = eye(4);
end
