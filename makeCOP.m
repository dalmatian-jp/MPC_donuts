function cop = makeCOP(x, u, params, pvstate)
    q1 = x(1);
    q2 = x(2);
    q1dot = x(3);
    q2dot = x(4);
    
    % Retrieve physical parameters from pvstate
    l1 = pvstate(1);    % Length of link 1 (m)
    s1 = pvstate(2);    % Distance from joint 1 to center of mass of link 1 (m)
    m1 = pvstate(3);    % Mass of link 1 (kg)
    J1 = pvstate(4);    % Moment of inertia of link 1 (kg*m^2)
    l2 = pvstate(5);    % Length of link 2 (m)
    s2 = pvstate(6);    % Distance from joint 2 to center of mass of link 2 (m)
    m2 = pvstate(7);    % Mass of link 2 (kg)
    J2 = pvstate(8);    % Moment of inertia of link 2 (kg*m^2)
    g  = pvstate(9);    % Acceleration due to gravity (m/s^2)
    q0 = 4 * pi / 180;
    l0 = 0.27;
    k1 = 0.468;
    k2 = 0.73;
    
    dxdt = doubleLinkDynamics(0, x, u, params);
    ddq1 = dxdt(3);
    ddq2 = dxdt(4);
    
    ut = (-(m1*k1^2 + m2*l1^2 + m2*k2^2)+2*m2*l1*k2*cos(q2)+(m1*l0*k1 + m2*l0*l1)*cos(q0)*cos(q1-q0)-m2*l0*k2*cos(q0)*cos(q2+q1-q0))*ddq1 ...
        +(-m2*k2^2 + m2*l1*k2*cos(q2)-m2*l0*k2*cos(q0)*cos(q2+q1-q0))*ddq2;
    
    fv = (m1*k1 +m2*l1)*cos(q1-q0)*ddq1...
        -(m2*k2 + m2*l2)*cos(q2 +q1 -q0)*(ddq2 + ddq1) +(m1+m2)*g;
    
    cop = ut/fv;
    end