function [Xpo, Ppo] = update_state(Ppr, Xpr, H, error, std_meas)
% Calculate Kalman Gain for the incoming sensor measurement
Kk = Ppr*H'*(H*Ppr*H' + std_meas^2)^(-1);

% Update estimate and covariance posteriors
Xpo = Xpr' + Kk*error;
Ppo = (eye(9) - Kk*H)*Ppr;

end