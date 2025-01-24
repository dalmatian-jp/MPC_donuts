function ekfDoubleLink()
    clc; close all;

    %% 1) Define system parameters
    systemParams = struct(...
        'm1', 25.00, ...
        'm2', 50.00, ...
        'L1', 0.90, ...
        'L2', 0.85, ...
        'I1', 1.00, ...
        'I2', 1.20, ...
        'com1', 0.45, ...
        'com2', 0.42, ...
        'g',  9.81);

    dt = 0.01;       % sampling time
    N  = 10;         % MPC horizon
    numSteps = 400;  % total simulation steps

    %% 2) Define Q, R for the MPC
    Q_mpc = diag([6457.5, 7801.7, 7838.3, 8860.4]);  % 4x4
    R_mpc = diag([0.02, 0.01]);                      % 2x2

    %% 3) Define Noise Covariances for EKF
    % Small process noise
    Q_ekf = 1e-6 * eye(4);
    % Measurement noise on all 4 states
    R_ekf = diag([1e-5, 1e-5, 1e-4, 1e-4]);

    %% 4) Initial True State + EKF estimate
    x_true = [0.0873; 0; 0; 0];         % slight forward lean
    x_hat  = x_true + 0.01*randn(4,1); % initial guess
    P      = 1e-3 * eye(4);           % initial covariance

    x_ref  = [0;0;0;0];
    
    %% 5) Data storage
    stateHistory     = zeros(numSteps+1,4);   % true
    estStateHistory  = zeros(numSteps+1,4);   % estimate
    stateHistory(1,:)= x_true';
    estStateHistory(1,:)= x_hat';

    %% 6) Main simulation loop
    for k=1:numSteps
        % ---- (A) Generate a measurement from the true state (with noise) ----
        z_meas = x_true + mvnrnd(zeros(1,4), R_ekf)';  % measure all 4 states

        % ---- (B) EKF Prediction step ----
        [x_pred, P_pred] = ekfPredict(x_hat, P, ...
                                      @(xx,uu) fEuler(xx,uu,systemParams,dt), ...
                                      @(xx,uu) A_jacobian(xx,uu,systemParams,dt), ...
                                      zeros(4,1), Q_ekf, ...  % no input noise
                                      mpcController(x_hat, x_ref, systemParams, N, Q_mpc, R_mpc, dt));
        % Here we used mpcController(...) as 'u', so we "pretend" to know the same control.

        % ---- (C) EKF Update step ----
        [x_hat_new, P_new] = ekfUpdate(x_pred, P_pred, ...
                                       z_meas, ...
                                       @(xx) h_meas(xx), ...
                                       @(xx) H_jacobian(xx), ...
                                       R_ekf);

        x_hat = x_hat_new;
        P     = P_new;

        % ---- (D) Use x_hat in MPC to get next control input ----
        u_opt = mpcController(x_hat, x_ref, systemParams, N, Q_mpc, R_mpc, dt);

        % ---- (E) Propagate the true system with real dynamics + process noise
        % We'll do a single Euler step for simplicity
        dxdt_true = doubleLinkDynamics(0, x_true, u_opt, systemParams);
        x_true    = x_true + dt*dxdt_true + mvnrnd(zeros(1,4), Q_ekf)';

        % ---- (F) Log data
        stateHistory(k+1,:)    = x_true';
        estStateHistory(k+1,:) = x_hat';
    end

    %% A) Angles and Velocities vs. Time
    timeVector = 0:numSteps;  % discrete time steps (if dt=0.01, actual time = timeVector*dt)
    
    figure('Name','Angle and Velocity','Color','white');
    
    % 1) Angles subplot
    subplot(2,1,1); 
    plot(timeVector, stateHistory(:,1), 'b-',  'LineWidth',2, 'DisplayName','q1 True'); hold on;
    plot(timeVector, estStateHistory(:,1),'b--','LineWidth',1.5,'DisplayName','q1 Est');
    plot(timeVector, stateHistory(:,2), 'r-',  'LineWidth',2, 'DisplayName','q2 True');
    plot(timeVector, estStateHistory(:,2),'r--','LineWidth',1.5,'DisplayName','q2 Est');
    xlabel('Time Step'); ylabel('Angle (rad)');
    title('q1, q2 (True vs. Est)');
    legend('Location','best'); grid on;
    
    % 2) Velocities subplot
    subplot(2,1,2);
    plot(timeVector, stateHistory(:,3), 'm-',  'LineWidth',2, 'DisplayName','dq1 True'); hold on;
    plot(timeVector, estStateHistory(:,3),'m--','LineWidth',1.5,'DisplayName','dq1 Est');
    plot(timeVector, stateHistory(:,4), 'g-',  'LineWidth',2, 'DisplayName','dq2 True');
    plot(timeVector, estStateHistory(:,4),'g--','LineWidth',1.5,'DisplayName','dq2 Est');
    xlabel('Time Step'); ylabel('Angular Velocity (rad/s)');
    title('dq1, dq2 (True vs. Est)');
    legend('Location','best'); grid on;

end

%% ------------------------------------------------------------------------
%% Helper function: discrete process update (Euler)
function x_next = fEuler(x, u, params, dt)
    dxdt = doubleLinkDynamics(0, x, u, params);
    x_next = x + dt*dxdt;
end

%% ------------------------------------------------------------------------
%% Jacobian of f w.r.t. x (for the prediction step)
function A = A_jacobian(x, u, params, dt)
    % Differentiate doubleLinkDynamics, then approximate discrete by (I + dt*A_cont).
    syms q1 q2 dq1 dq2 tau1 tau2 real
    state = [q1; q2; dq1; dq2];
    ctrl  = [tau1; tau2];

    dxdt_sym = doubleLinkDynamics(0, state, ctrl, params);
    A_sym    = jacobian(dxdt_sym, state);

    A_cont = double(subs(A_sym, [state; ctrl], [x; u]));
    A      = eye(4) + dt*A_cont;
end

%% ------------------------------------------------------------------------
%% Measurement function h(x)
function z = h_meas(x)
    % We assume we measure the entire state directly
    z = x;
end

%% ------------------------------------------------------------------------
%% Jacobian of h w.r.t. x
function H = H_jacobian(x)
    % If z = x, then H = Identity(4x4)
    H = eye(4);
end
