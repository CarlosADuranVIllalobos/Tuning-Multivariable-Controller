% To simulate the Fed-batch mammalian cell biprocess 
% Copyright Carlos Alberto Duran Villalobos August 2021- University of Manchester.
% Last modified: March , 2021
% See Craven, Whelan, Glennon 2014 Paper - Journal of process control 24, 344-357

function [X,unom_Fs,yest] = mammalianFBLoopMPC(Y0,tt,step,us,contr_param)
% INPUTS
% Y0 - initial conditions 
% tt - total simulation time
% step - sampling time (timestep)
% us - trajectories to follow 
% OUTPUTS
% X - matrix with batch data at each timestep

t=(0:step:tt)'; % get simulation time steps up to end time
cp = [contr_param.ustart/step:contr_param.step/step:contr_param.bt/step]; %Controlpoints
mammalianFBParams; % call script to get model parameters

% get initial conditions from Y0
Xt = Y0(1); % total cell density
Xv = Y0(2); % viable cells
Xd = Y0(3); % dead cells
G = Y0(4);  % glucose concentration
Q = Y0(5);  % glutamie concentration 
L = Y0(6);  % lactate concentration
A = Y0(7);  % ammonia concentration
Ci = Y0(8); % unknown inhibitor concentration
V = Y0(9);  % volume
F = us(1,:); % glucose Flowrate set point
% calculate initial u and kd (growth rate and death rate)
u = (umax*G*Q*(Cstari-Ci)) / (Cstari*(Kg+G)*(Kq+Q)*((L/Kl)+1)*((A/Ka)+1)); % growth rate
Kd = Kdmax*(Ku /(u+Ku)); % death rate

step1 = 0.01; % step size for ODE solver
X = [];
for ts = 1:1:(tt/step)
    %V = F(ts) + V;
    if any(ts==cp)% Add MPC controller
        [Recipe_Fs,yest] = fast_controller(X,contr_param,ts*step,F);
        F = Recipe_Fs;
    end
    
    ODEinp = [Xt Xv Xd G Q L A Ci V];
    [T1,Y1] = ode45('mammalianFBODEs',t(ts):step1:t(ts+1),ODEinp,[],[F(ts) V step1]);
    [m1,n1] = size(Y1);
    Y1 = Y1(m1,:);
    X = [X; Y1 F(ts)];
    Xt = Y1(1);
    Xv = Y1(2);
    Xd = Y1(3);
    G = Y1(4);
    Q = Y1(5);
    L = Y1(6);
    A = Y1(7);
    Ci = Y1(8);
    V = Y1(9);
    %F = Y1(10);    
end
unom_Fs=F;   
    
    
