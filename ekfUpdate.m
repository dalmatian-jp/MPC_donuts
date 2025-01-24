function [xPost, PPost] = ekfUpdate(xPred, PPred, z, hFun, HJacFun, R)
% ekfUpdate: One-step EKF measurement update
%   xPred, PPred: predicted state and covariance
%   z           : actual measurement
%   hFun(x)     : measurement function
%   HJacFun(x)  : returns H = dh/dx
%   R           : measurement noise cov

    % 1) Compute expected measurement
    zPred = hFun(xPred);

    % 2) Linearize measurement model
    H = HJacFun(xPred);

    % 3) Kalman Gain
    S = H*PPred*H' + R;
    K = PPred*H'*(S \ eye(size(S)));  % or inv(S) * for smaller dims

    % 4) Update
    xPost = xPred + K*(z - zPred);
    PPost = (eye(size(PPred)) - K*H)*PPred;
end