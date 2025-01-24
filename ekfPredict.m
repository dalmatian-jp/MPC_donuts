function [xPred, PPred] = ekfPredict(xPrev, PPrev, fFun, AJacFun, wMean, Q, u)
% ekfPredict: One-step EKF prediction
%   xPrev, PPrev: previous state estimate and cov
%   fFun(x,u)   : discrete-time process model
%   AJacFun(x,u): returns A = df/dx at x
%   wMean       : mean of process noise (often zero)
%   Q           : process noise cov
%   u           : known control input

    % 1) Predict state
    xPred = fFun(xPrev, u) + wMean;  % e.g. Euler step of the dynamics

    % 2) Predict covariance
    A     = AJacFun(xPrev, u);
    PPred = A*PPrev*A' + Q;
end

