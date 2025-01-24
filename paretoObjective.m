function objectives = paretoObjective(params, x0, x_ref, N, dt, systemParams)
    % Extract Q and R from 'params'
    Q = diag(params(1:4));  % e.g. q1..q4
    R = diag(params(5:6));  % r1..r2

    numSteps = 50;
    x        = x0;
    totalErr = 0;
    totalU   = 0;

    for k = 1:numSteps
        try
            u_opt = mpcController(x, x_ref, systemParams, N, Q, R, dt);
        catch ME
            warning("mpcController error at step %d: %s", k, ME.message);
            u_opt = [0;0];
        end

        % Integrate one step
        [~, X_sim] = ode45(@(t,xx) doubleLinkDynamics(t, xx, u_opt, systemParams), [0 dt], x);
        x = X_sim(end,:)';

        % Accumulate
        totalErr = totalErr + norm(x - x_ref)^2;
        totalU   = totalU + norm(u_opt)^2;
    end

    rmse   = sqrt(totalErr / numSteps);
    energy = totalU   / numSteps;
    
    % Debugging objectives
    disp(['RMSE: ', num2str(rmse), ', Energy: ', num2str(energy)]);
    
    objectives = [rmse, energy];
end


