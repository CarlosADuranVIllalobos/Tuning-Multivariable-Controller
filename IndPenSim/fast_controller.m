
function [Recipe_Fs_sp,Recipe_PAA_sp] = fast_controller(Xref,c_param,cp,unom_Fs,unom_PAA)
% Function implements MPC control inside an ODE loop
% Copyright Carlos Alberto Duran Villalobos August 2021 - University of Manchester.

Xref.pH.y = -log(Xref.pH.y) / log(10); % convert to pH from H+ concentration
Xref.Q.y = Xref.Q.y / 1000; % convert heat from Qrxn to kcal

%unom_Fs=c_param.Recipe_Fs_sp;
%unom_PAA=c_param.Recipe_PAA_sp;


[xc,yc]=convfactor(Xref,c_param.ustart); %convert values for PLS
xc=[c_param.x0 xc];
%orgplsmpc for M constant, orgplsmpcM for M tuned
[Xtrain,xp,Porg,Worg,Eorg,xmean2,xstd2]=orgplsmpcM(c_param.X2,c_param.P,c_param.W,c_param.E,c_param.xmean,...
    c_param.xstd,c_param.num_x,c_param.num_u,c_param.bt,c_param.ustart,xc,cp,c_param.invar);  %organize variables for MPC
xk=[xp unom_Fs(cp:end) unom_PAA(cp:end) ]; %known variables
np=length(xp);
nk=length(xk);
nu=length(unom_Fs(cp:end))+length(unom_PAA(cp:end));
xk = (xk-xmean2(1:nk))./xstd2(1:nk);

%Feed constraints
lower_Fs=5*ones(1,c_param.timeinterval-cp+1); %Lower limit Fs
upper_Fs=180*ones(1,c_param.timeinterval-cp+1); %Upper limit Fs
lower_PAA=0.01*ones(1,c_param.timeinterval-cp+1); %Lower limit PAA
upper_PAA=15*ones(1,c_param.timeinterval-cp+1); %Upper limit PAA
lower=[lower_Fs lower_PAA];
upper=[upper_Fs upper_PAA];
%% Change here for M!!!!!
%Begin optimisation
M=diag(c_param.M(cp*c_param.num_u-1:end));%MVT control tunning

vol_con=29000; %constraint on volume
ysppls=(c_param.ysp-c_param.ymean)./c_param.ystd;%
[ufuture,H,Fu,tk]=mvtoptujtje(Xtrain,ysppls,xk,np,nu,nk,Eorg,c_param,Porg,Worg,xmean2,xstd2,upper,lower,M,vol_con);

%Run batch simulation with the optimised MVT
unom_Fs(cp:end)=ufuture(1:c_param.timeinterval-cp+1)';
unom_PAA(cp:end)= ufuture(c_param.timeinterval-cp+2:end)';

%% Continue
Recipe_Fs_sp = unom_Fs;
Recipe_PAA_sp= unom_PAA;
end


