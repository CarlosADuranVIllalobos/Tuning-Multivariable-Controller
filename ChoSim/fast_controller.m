
function [Recipe_Fs_sp] = fast_controller(Xref,c_param,cp,unom_Fs)
% Function implements MPC control inside an ODE loop- only runs 1 simulation per run
% For indpensim_run.m from Stephen Goldrick September 2014
% Copyright Carlos Alberto Duran Villalobos August 2021 - University of Manchester.

[xc,yc]=convfactor(Xref,c_param.ustart); %convert values for PLS
xc=[c_param.x0(4) c_param.x0(5) c_param.x0(6) c_param.x0(7) c_param.x0(9) xc];
%orgplsmpc for M constant, orgplsmpcM for M tuned
[Xtrain,xp,Porg,Worg,Eorg,xmean2,xstd2]=orgplsmpcM(c_param.X2,c_param.P,c_param.W,c_param.E,c_param.xmean,...
    c_param.xstd,c_param.num_x,c_param.num_u,c_param.bt,c_param.ustart,xc,cp,c_param.invar);  %organize variables for MPC
xk=[xp unom_Fs(cp:end) ]; %known variables
np=length(xp);
nk=length(xk);
nu=length(unom_Fs(cp:end));
xk = (xk-xmean2(1:nk))./xstd2(1:nk);

%Feed constraints
lower=0*ones(1,c_param.timeinterval-cp+1); %Lower limit Fs
upper=0.04*ones(1,c_param.timeinterval-cp+1); %Upper limit Fs

%Control tuning matrix
M=diag(c_param.M(cp*c_param.num_u:end));%MVT control tunning

vol_con =15-c_param.x0(9); %constraint on volume
ysppls=(c_param.ysp-c_param.ymean)./c_param.ystd;%
[ufuture,H,Fu,tk]=mvtoptujtje(Xtrain,ysppls,xk,np,nu,nk,Eorg,c_param,Porg,Worg,xmean2,xstd2,upper,lower,M,vol_con);
%Run batch simulation with the optimised MVT
unom_Fs(cp:end)=ufuture(1:c_param.timeinterval-cp+1)';

%% Continue
Recipe_Fs_sp = unom_Fs;
end


