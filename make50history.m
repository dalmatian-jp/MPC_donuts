clc; close all;

function  [estSmall,s_StateHistory,s_uHistory] = for50(numSteps)
    % twoPerturbationsEKF:
    %   Demonstrates running the double-link EKF+MPC simulation for two
    %   different initial conditions ("small" vs "large" perturbation),
    %   then plots only the ESTIMATED phase portraits for comparison,
    %   with a progress bar for each run.
    
        clc; close all;
    
        %% 1) Define system parameters
        systemParams = struct(...
            'm1',   28.00, ...
            'm2',   53.00, ...
            'L1',   0.90, ...
            'L2',   0.88, ...
            'I1',   9.21, ...
            'I2',   5.35, ...
            'com1', 0.58, ...
            'com2', 0.32, ...
            'g',    9.81);
    
        %% 2) MPC parameters
        dt       = 0.01;                % sampling time
        N        = 10;                  % MPC horizon
        Q_mpc    = diag([6457.5, 7801.7, 7838.3, 8860.4]); % 4x4
        R_mpc    = diag([0.02, 0.01]);  % 2x2
        % numSteps = 100;                % how many discrete steps to simulate
    
        %% 3) EKF noise covariances
        Q_ekf = 1e-6 * eye(4);
        R_ekf = diag([1e-4, 1e-4, 1e-4, 1e-4]);
    
        %% 4) Define two initial states (perturbations)
        % "Small" perturbation
        x0_small = [0.0873; 0; 0; 0];  % from the research article
        % "Large" perturbation
        % x0_large = [0.20; 0; 0; 0];    % bigger lean (0.2618 in research article)
    
        x_ref = [0; 0; 0; 0];          % same reference for both
    
        %% 5) Run simulation for each perturbation (with progress bars)
        [estSmall,s_StateHistory,s_uHistory]  = runEKFsimulation(...
            x0_small, x_ref, systemParams, dt, N, Q_mpc, R_mpc, Q_ekf, R_ekf, numSteps, ...
            'Small Perturbation'  ... % just to label the waitbar
        );
        % [estLarge,l_StateHistory,l_uHistory] = runEKFsimulation(...
        %     x0_large, x_ref, systemParams, dt, N, Q_mpc, R_mpc, Q_ekf, R_ekf, numSteps, ...
        %     'Large Perturbation'  ...
        % );
    
        %% 6) Plot a single-phase portrait for ANKLE (q1,dq1) with 2 trajectories
        figure('Name','Ankle Phase Portrait','Color','white');
        plot(estSmall(:,1), estSmall(:,3), 'b-', 'LineWidth',2, ...
            'DisplayName','Small Perturbation'); hold on;
        % plot(estLarge(:,1), estLarge(:,3), 'r--','LineWidth',2, ...
        %     'DisplayName','Large Perturbation');
        % xlabel('q_1 (rad)'); ylabel('dq_1 (rad/s)');
        % title('Estimated Ankle Phase Portrait');
        % legend('Location','best'); grid on;
    
        %% 7) Plot a single-phase portrait for HIP (q2,dq2) with 2 trajectories
        figure('Name','Hip Phase Portrait','Color','white');
        plot(estSmall(:,2), estSmall(:,4), 'b-', 'LineWidth',2, ...
            'DisplayName','Small Perturbation'); hold on;
        % plot(estLarge(:,2), estLarge(:,4), 'r--','LineWidth',2, ...
        %     'DisplayName','Large Perturbation');
        % xlabel('q_2 (rad)'); ylabel('dq_2 (rad/s)');
        % title('Estimated Hip Phase Portrait');
        % legend('Location','best'); grid on;
    
    end  % end of twoPerturbationsEKF
    
    
    %% ========================================================================
    function [estStateHistory,StateHistory,uHistory] = runEKFsimulation(...
                            x0, x_ref, systemParams, dt, N, Q_mpc, R_mpc, ...
                            Q_ekf, R_ekf, numSteps, runLabel)
    % runEKFsimulation:
    %   Runs an EKF+MPC loop for `numSteps` starting from initial state x0.
    %   Only returns the ESTIMATED states. 
    %   A waitbar (progress bar) is displayed during the simulation.
    %
    %   Input 'runLabel' is just used for the waitbar title.
    
        % For demonstration, we create a "true" state to propagate with noise.
        x_true = x0;
    
        % EKF state/cov
        x_hat = x_true + 0.01*randn(4,1);   % slightly perturbed initial guess
        P     = 1e-3*eye(4);               % initial cov
    
        % Preallocate log
        estStateHistory = zeros(numSteps+1,4);
        estStateHistory(1,:) = x_hat';
        StateHistory = zeros(numSteps+1,4);
        StateHistory(1,:) = x_true';
        uHistory = zeros(numSteps+1,2);
    
        % Create a waitbar for progress
        wb = waitbar(0, ['Running EKF Simulation: ' runLabel]);
    
        for k = 1:numSteps
    
            % 1) Simulate a measurement from the true state
            z_meas = x_true + mvnrnd(zeros(1,4), R_ekf)';  % measure all 4 states
    
            % 2) EKF Predict: use the same control law to get 'u_pred'
            u_pred = mpcController(x_hat, x_ref, systemParams, N, Q_mpc, R_mpc, dt);
            [xPred, PPred] = ekfPredict(x_hat, P, ...
                @(xx,uu) fEuler(xx,uu,systemParams,dt), ...
                @(xx,uu) A_jacobian(xx,uu,systemParams,dt), ...
                zeros(4,1), Q_ekf, u_pred);
    
            % 3) EKF Update
            [x_hat, P] = ekfUpdate(xPred, PPred, ...
                z_meas, @(xx) h_meas(xx), @(xx) H_jacobian(xx), R_ekf);
    
            % 4) Use x_hat for the actual control
            u_opt = mpcController(x_hat, x_ref, systemParams, N, Q_mpc, R_mpc, dt);
    
            % 5) Propagate the TRUE system with noise
            % dxdt_true = doubleLinkDynamics(0, x_true, u_opt, systemParams);
            dxdt_true = doubleLinkDynamics_takami(0, x_true, u_opt, systemParams);
            % Simple Euler step + process noise
            x_true = x_true + dt*dxdt_true + mvnrnd(zeros(1,4), Q_ekf)';
    
            % 6) Log the estimated state
            estStateHistory(k+1,:) = x_hat';
            StateHistory(k+1,:) = x_true';
            uHistory(k+1,:) = u_opt';
    
            % 7) Update progress bar
            waitbar(k/numSteps, wb);
        end
    
        % Close the progress bar
        close(wb);
    
    end
    
    
    %% ========================================================================
    % The next few subfunctions are for the EKF approach.
    
    function [xPred, PPred] = ekfPredict(xPrev, PPrev, fFun, AJacFun, wMean, Q, u)
    % EKF Prediction step
        xPred = fFun(xPrev, u) + wMean;
        A     = AJacFun(xPrev, u);
        PPred = A*PPrev*A' + Q;
    end
    
    function [xPost, PPost] = ekfUpdate(xPred, PPred, z, hFun, HJacFun, R)
    % EKF Update step
        zPred = hFun(xPred);
        H     = HJacFun(xPred);
        S     = H*PPred*H' + R;
        K     = PPred*H'*(S \ eye(size(S)));
        xPost = xPred + K*(z - zPred);
        PPost = (eye(size(PPred)) - K*H)*PPred;
    end
    
    %% ========================================================================
    % The "discrete" process model = Euler step of doubleLinkDynamics
    function x_next = fEuler(x, u, params, dt)
        % dxdt  = doubleLinkDynamics(0, x, u, params);
        dxdt  = doubleLinkDynamics_takami(0, x, u, params);

        x_next= x + dt*dxdt;
    end
    
    % Jacobian of the discrete model wrt x, using symbolic or numerical partials
    function A = A_jacobian(x, u, params, dt)
        syms q1 q2 dq1 dq2 tau1 tau2 real
        state = [q1; q2; dq1; dq2];
        ctrl  = [tau1; tau2];
    
        % dxdt_sym = doubleLinkDynamics(0, state, ctrl, params);
        dxdt_sym = doubleLinkDynamics_takami(0, state, ctrl, params);
        A_sym    = jacobian(dxdt_sym, state);
    
        A_cont = double(subs(A_sym, [state; ctrl], [x; u]));
        A      = eye(4) + dt*A_cont;
    end
    
    % Measurement function: We measure the entire 4D state with noise
    function z = h_meas(x)
        z = x;
    end
    
    % Jacobian of measurement function h wrt x
    function H = H_jacobian(~)
        H = eye(4);
    end
    

for i = 101:150
    [estSmall, s_StateHistory, s_uHistory] = for50(3000);
    allDataMatrix = [estSmall, s_StateHistory, s_uHistory];
                    %  estLarge, l_StateHistory, l_uHistory];

    % CSVファイルに書き出し
    filename = sprintf('new50\\EKF_Simulation_AllData_%02d.csv', i);
    writematrix(allDataMatrix, filename);

    disp(['シミュレーションデータを ' filename ' に保存しました。']);
end
