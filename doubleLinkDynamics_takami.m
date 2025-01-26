function dxdt = doubleLinkDynamics_takami(t, x, u, params)
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

    r1 = 0.64*L1;
    r2 = 0.36*L2;

    % Unpack states
    q1 = x(1); q2 = x(2); % Angles
    dq1 = x(3); dq2 = x(4); % Angular velocities

    q1_new = q1 + pi/2;

    a = I1 + I2 + m1*r1^2 + m2*(L1^2 + r2^2);
    b = m2*L1*r2;
    d = I2 + m2*r2^2;
    
    s1 = sin(q1_new);
    s2 = sin(q2);
    c1 = cos(q1_new);
    c2 = cos(q2);
    c12 = cos(q1_new + q2);

    % Compute the Mass Matrix (M)
    M11 = a+2*b*c2;
    % M12
    M12 = d+b*c2;

    % M行列
    M = [ M11, M12;
          M12, d];


    % Compute the Coriolis Matrix (C)
    C = [-b*s2*dq2, -b*s2*(dq1+dq2);
           b*s2*dq2, 0];

    % Compute the Gravity Vector (F)
    G1 = -g*((m1*r1 + m2*L1)*c1 + m2*r2*c12);
    G2 = -g*m2*r2*c12;
    G = [G1; G2];

    % Solve for angular accelerations
    ddq = (M + 1e-6 * eye(2)) \ (u - C*[dq1;dq2] - G);

    % Return state derivatives
    dxdt = [dq1; dq2; ddq];
end
