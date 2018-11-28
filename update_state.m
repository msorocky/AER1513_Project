function [Xpo, Ppo] = update_state(Ppr, Xpr, H, error, R)
% Calculate Kalman Gain for the incoming sensor measurement
Kk = Ppr*H'*(H*Ppr*H' + R)^(-1);

% Update estimate and covariance posteriors
Xpo = Xpr' + Kk*error;
Ppo = (eye(9) - Kk*H)*Ppr;

% Enforce symmetry of covariance matrix
Ppo = 0.5*(Ppo+Ppo');

end