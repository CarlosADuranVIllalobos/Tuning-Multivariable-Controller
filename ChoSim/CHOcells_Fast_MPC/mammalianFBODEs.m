% Provide ODEs for mammalian fed batch - see mammalianFBLoop
% Copyright Carlos Alberto Duran Villalobos August 2021- University of Manchester.
% Last modified: March , 2021

function dy = mammalianFBODEs(t,y,FLAG,inp);
% Provides ODEs for mammalian Fed batch loop
% INPUTS:
% t - time
% y - initial conditions for ODE states
% inp - inputs (F,V)
% OUTPUTS:
% dy - change in y to be calculated by solver

% call in model parameters from script
mammalianFBParams;
F = inp(1); % get flowrate and volume from input
V = inp(2);
step1 = inp(3);
u = (umax*y(4)*y(5)*(Cstari-y(8))) / (Cstari*(Kg+y(4))*(Kq+y(5))*((y(6)/Kl)+1)*((y(7)/Ka)+1));
Kd = Kdmax*(Ku /(u+Ku));
% get change in ys
dy(1) = u*y(2)-Klysis*y(3) - y(1)*(F/V); % total cell density
dy(2) = (u-Kd)*y(2) - y(2)*(F/V); % viable cells
dy(3) = Kd*y(2) - Klysis*y(3) - y(2)*(F/V); % cell death
dy(4) = (F/V)*(Sg-y(4)) + ((-u/Yxg)-mg)*y(2); % change in glucose
%dy(4) = ((-u/Yxg)-mg)*y(2);
dy(5) = (F/V)*(Sq-y(5)) + ((-u/Yxq)-mq)*y(2) - Kd*y(5); % glutamine change
%dy(5) = ((-u/Yxg)-mq)*y(2) - Kdq*y(5); % glutamine change
dy(6) = -Ylg*((-u/Yxg)-mg)*y(2) - (F/V)*y(6); % lactate change 
dy(7) = -Yaq*((-u/Yxq)-mq)*y(2) + Kdq*y(5) - (F/V)*y(7); % ammonia change
dy(8) = qi*y(2) - (F/V)*y(8); % inhibitor concentration change
dy(9) = F;
dy=[dy(1);dy(2);dy(3);dy(4);dy(5);dy(6);dy(7);dy(8);dy(9)];

end

