function [u] = PIDSimple3(uk1, ek, ek1, yk, yk1, yk2, u_min, u_max, Kp, Ti, Td, h)
%% PIDsimple2.m
% Simple PID controller in incremental form. The derivative component is
% applied only to the output. Controller also checks for control signal
% saturation. Integral and derivative time constants must be defined in the
% same time unit as the sampling period.
% 
% [u] = PIDSimple(uk1, ek, ek1, yk, yk1, yk2, u_min, u_max, Kp, Ti, Td, h)
%
% Parameters:
%   uk1 - previous control input
%   ek - current error
%   ek1 - previous error
%   yk - current PV
%   yk1 - previous PV
%   yk2 - previous previous PV
%   Kp - proportional gain
%   Ti - integral time constant
%   Td - derivative time constant
%   h - sampling time


%%
% 
% $$u_{k} = u_{k-1} + K_p \left[ e_k - e_{k-1} + \frac{e_k}{T_i} + T_d(y_k-2y_{k-1}+y_{k-2})\right]$$
% 
% proportional component
P = ek - ek1;

% checks if the integral time constant is defined
if(Ti>1e-7)
    I = ek*h/Ti;
else
    I = 0;
end

% derivative component
if(Td > 0.001)
    D = -Td/h * (yk - 2*yk1 + yk2);
else
    D = 0;
end


% computes and saturates the control signal
u = uk1 + Kp*(P + I + D);
if(u>u_max)
    u = u_max;
end
if(u<u_min)
    u = u_min;
end

end
%% Copyright
% Stephen Goldrick and Andrei Stefan
% University of Manchester, Newcastle University and Perceptive Engineering
%
% All rights reserved. Copyright (c) Newcastle University, University of Manchester and Perceptive Engineering.
