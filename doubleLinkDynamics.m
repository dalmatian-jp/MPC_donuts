function dxdt = doubleLinkDynamics(t, x, u, params)
    %#codegen
    % doubleLinkDynamics: Computes the state derivatives for a double inverted pendulum.
    %
    % Inputs:
    %   t      - Time (not used in dynamics directly but kept for compatibility with ODE solvers).
    %   x      - State vector [q1; q2; dq1; dq2].
    %   u      - Control input vector [tau1; tau2].
    %   params - Structure containing system parameters:
    %            params.m1, params.m2 - Masses of the links.
    %            params.L1, params.L2 - Lengths of the links.
    %            params.I1, params.I2 - Moments of inertia.
    %            params.g             - Gravitational acceleration.
    %
    % Outputs:
    %   dxdt   - State derivatives [dq1; dq2; ddq1; ddq2].

    % Unpack parameters
    m1 = params.m1; m2 = params.m2;
    L1 = params.L1; L2 = params.L2;
    I1 = params.I1; I2 = params.I2;
    g = params.g;

    % Unpack states
    q1 = x(1); q2 = x(2); % Angles
    dq1 = x(3); dq2 = x(4); % Angular velocities

    % Compute the Mass Matrix (M)
    M11 = m1*L1^2/3 + m2*(L1^2 + (L2^2)/3 + L1*L2*cos(q2)) + I1 + I2;
    M12 = m2*((L2^2)/3 + L1*L2*cos(q2)) + I2;
    M21 = M12;
    M22 = m2*(L2^2)/3 + I2;
    M = [M11, M12; M21, M22];

    % Compute the Coriolis Matrix (C)
    C1 = -m2*L1*L2*sin(q2)*dq2^2 - 2*m2*L1*L2*sin(q2)*dq1*dq2;
    C2 = m2*L1*L2*sin(q2)*dq1^2;
    C = [C1; C2];

    % Compute the Gravity Vector (F)
    F1 = (m1*L1/2 + m2*L1)*g*sin(q1) + m2*(L2/2)*g*sin(q1 + q2);
    F2 = m2*(L2/2)*g*sin(q1 + q2);
    F = [F1; F2];

    % Solve for angular accelerations
    ddq = (M + 1e-6 * eye(2)) \ (u - C - F);

    % Return state derivatives
    dxdt = [dq1; dq2; ddq];
end
