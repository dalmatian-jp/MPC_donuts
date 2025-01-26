function u_opt = mpcController(x, x_ref, params, N, Q, R, dt)
    %#codegen
    % mpcController: Computes the optimal control inputs using MPC.

    % System dimensions
    nx = length(x);       % Number of states
    nu = 2;              % Number of control inputs (torques)
    
    % Check dimension
    if nx ~= 4
        error('State dimension must be 4 for double-link system.');
    end

    % Linearize system (local A, B matrices)
    [A, B] = linearizeSystem(x, params, dt);

    % Build prediction matrices
    [Phi, Gamma] = buildPredictionMatrices(A, B, N);

    % Debugging prediction matrices
    % disp('Prediction matrix dimensions:');
    % disp(['Phi size: ', mat2str(size(Phi))]);
    % disp(['Gamma size: ', mat2str(size(Gamma))]);

    % Build cost function matrices
    H = blkdiag(kron(eye(N), Q), kron(eye(N), R)); % Quadratic cost
    f = [repmat(-Q * x_ref, N, 1); zeros(N * nu, 1)]; % Penalize deviation from reference

    % Debugging cost function matrices
    % disp('Cost function dimensions:');
    % disp(['H size: ', mat2str(size(H))]);
    % disp(['f size: ', mat2str(size(f))]);

    % Equality constraints (system dynamics)
    [Aeq, beq] = buildEqualityConstraints(Phi, Gamma, x, nx, nu, N);

    % Debugging equality constraints
    % disp('Equality constraints dimensions:');
    % disp(['Aeq size: ', mat2str(size(Aeq))]);
    % disp(['beq size: ', mat2str(size(beq))]);

    % Bounds for decision variables
    angle_min = [-0.35; -0.53]; % Minimum joint angles (rad)
    angle_max = [0.53; 0.87];   % Maximum joint angles (rad)

    % Suppose you also want to bound the velocities:
    vel_min   = [-0.5; -0.5];         % [dq1_min; dq2_min] example
    vel_max   = [ 0.5;  0.5];

    torque_min = [-20; -40];    % Minimum torques (Nm)
    torque_max = [20; 40];      % Maximum torques (Nm)
    
    % For each state vector step, we have 4 elements: [q1,q2,dq1,dq2].
    lb_x = [angle_min(1); angle_min(2); vel_min(1); vel_min(2)];
    ub_x = [angle_max(1); angle_max(2); vel_max(1); vel_max(2)];

    % Replicate for N steps
    lb_states = repmat(lb_x, N, 1);
    ub_states = repmat(ub_x, N, 1);

    % For the inputs
    lb_controls = repmat(torque_min, N, 1);
    ub_controls = repmat(torque_max, N, 1);

    lb = [lb_states; lb_controls];
    ub = [ub_states; ub_controls];


    % Debugging bounds
    % disp('Bounds dimensions:');
    % disp(['lb size: ', mat2str(size(lb))]);
    % disp(['ub size: ', mat2str(size(ub))]);

    % Solve QP
    options = optimoptions('quadprog', 'Display', 'off');
    [z, ~, exitflag] = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], options);
    
    if exitflag ~= 1
        warning('QP infeasible or failed. Returning zero control.');
        u_opt = zeros(nu,1);
        return;
    end

    % Extract the first control input from the solution
    idx_u_start = N*nx + 1;    % states occupy the first N*nx
    idx_u_end   = N*nx + nu;   % first 2 controls
    u_opt = z(idx_u_start : idx_u_end);
    % disp(u_opt)
end

function [A, B] = linearizeSystem(x, params, dt)
    % Linearize the system around the current state using Jacobians
    % Inputs:
    %   x - Current state vector [q1; q2; dq1; dq2]
    %   params - System parameters
    %   dt - Sampling time
    % Outputs:
    %   A - State transition matrix
    %   B - Input matrix

    % Define symbolic variables for Jacobian computation
    syms q1 q2 dq1 dq2 tau1 tau2 real
    state = [q1; q2; dq1; dq2];
    ctrl = [tau1; tau2];

    dxdt_sym = doubleLinkDynamics(0, state, ctrl, params);
    % dxdt_sym = doubleLinkDynamics_takami(0, state, ctrl, params);
    A_sym = jacobian(dxdt_sym, state);
    B_sym = jacobian(dxdt_sym, ctrl);

    A_cont = double(subs(A_sym, [state; ctrl], [x; 0;0]));
    B_cont = double(subs(B_sym, [state; ctrl], [x; 0;0]));

    % Discrete approximation
    A = eye(4) + dt*A_cont;
    B = dt*B_cont;
end


function [Phi, Gamma] = buildPredictionMatrices(A, B, N)
    % Build state transition and input matrices over the prediction horizon
    [nx, nu] = size(B);

    Phi   = zeros(nx*N, nx);
    Gamma = zeros(nx*N, nu*N);

    for k = 1:N
        Phi((k-1)*nx + 1 : k*nx, :) = A^k;
        for j = 1:k
            Gamma((k-1)*nx + 1 : k*nx, (j-1)*nu + 1 : j*nu) = A^(k-j)*B;
        end
    end
end

function [Aeq, beq] = buildEqualityConstraints(Phi, Gamma, x0, nx, nu, N)
    % Build equality constraints for system dynamics

    % We want x(k) = A^k x0 + sum_{j=1 to k}(A^(k-j)*B*u(j)).
    % z = [ x(1); x(2); ... x(N); u(1); ... u(N)] 
    % => total length = (nx+nu)*N
    %
    % Aeq * z = beq enforces x(k) relationships for k=1..N.

    % Number of states stacked = nx*N
    % Number of inputs stacked = nu*N
    totalVars = (nx+nu)*N;

    Aeq = zeros(nx*N, totalVars);
    beq = zeros(nx*N, 1);

    % For each k in [1..N], the row block is (k-1)*nx+1 : k*nx.
    for k = 1:N
        row_ix = (k-1)*nx + (1:nx);
        % The sub-block from Phi, Gamma
        Phi_k   = Phi(row_ix, :);
        Gamma_k = Gamma(row_ix, :);

        % Indices for x(k) in z
        xk_start = (k-1)*nx + 1;
        xk_end   = k*nx;

        % Indices for the entire input block in z
        u_start = nx*N + 1;
        u_end   = nx*N + nu*N;

        % Place +I for x(k)
        Aeq(row_ix, xk_start:xk_end) = eye(nx);

        % Move x(k) to LHS => x(k) - [A^k x0 + sum(...) ] = 0
        beq(row_ix) = Phi_k * x0; % This is the known offset from x0

        % We must subtract the input contributions => -Gamma_k * U
        % But "Gamma_k" is  (nx x nu*N), we place it in Aeq with a negative sign
        Aeq(row_ix, u_start:u_end) = -Gamma_k;
    end
end