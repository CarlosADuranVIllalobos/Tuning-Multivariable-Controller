%Function to optimise an MVT
%Copyright Carlos Alberto Duran Villalobos July 2017-University of Manchester

function [Ufuture,H,Fu,yest,tk]=mvtoptu1(Xtrain,ysp,Xc,np,nu,nb,E,T,Porg,Worg,C,xmean2,xstd2,ymean,ystd,...
    upper,lower,M,vol_con)
% Xtrain: Training values for the PLS model 
% ysp: set-point for y
% Xc: nominal vector of predictors
% np: length of past predictors
% nu: length of predictors for the MVT
% nb: length of past predictors + MVT
% E: Residuals PLS
% T: Scores PLS
% Porg: X Loadings PLS organised by control points
% Worg: Weights PLS organised by control points
% C: Y LoadingsPLS  
% xmean2: mean of predictors organised by control points
% xstd2: std of predictors organised by control points
% ymean: mean of outputs
% ystd: std of outputs
% upper: Upper constraints MVT
% lower: lower constraints MVT
% M: matrix of weights for the optimisation problem
% vol_con: Volume constraint

%% Initialisation
X1 = Xc(1,1:np)';
Xp= Xc(1,1:np+nu)';
W1 = Worg(1:np,:);

Ppu = Porg(1:np+nu,:);
Pu = Porg(np+1:np+nu,:);
Pf = Porg(np+nu+1:end,:);

Wf = Worg(np+nu+1:end,:);
Wu = Worg(np+1:np+nu,:);
Wpu = Worg(1:np+nu,:);

U_future = Xc(:,np+1:np+nu)';

%Mising data by Trimed Score Regression
Xpu = Xtrain(:,1:np+nu);
Xf = Xtrain(:,np+nu+1:end);
beta=Wpu/(Wpu'*Xpu'*Xpu*Wpu)*Wpu'*Xpu'*Xtrain*Worg;
t_hat=Xp'*beta;
X2=(t_hat*Pf')';
D=beta*Pf'*Wf;
Du=D(np+1:np+nu,:)';
%Mising data by Known Data Regression
%       theta=T'*T;
%      % Beta=(Ppu*theta*Ppu')\Ppu*theta*Pf'*Wf+Wpu; %PLS
%       Beta=theta*Ppu'/(Ppu*theta*Ppu');
%       t_hat=Beta*Xp;
%Missing data by PMP
%       D =((Ppu'*Ppu)\Ppu');
%       Du=D(:,np+1:np+nu);
%     X2=Pf*D*Xp;

% Optimisation Routine tracking YSP with TSR
Xc=[Xp' X2'];
tp=X1'*W1;
tu=U_future'*Wu;
tf=X2'*Wf;

    function f = myfun(x,H,f)
        f= 0.5*x'*H*x + f'*x;
    end
    function [c,ceq] = mycon(x,param)
        % Calculate indices
        if (param.hard.Jt>0)
            T=param.tp+x'*param.td;
            Jt=T/param.sa*T'*param.Jt;
        end
        if (param.hard.Je>0)
            Ep=param.ep;
            Eu=x'*param.eu;
            Ef=x'*param.ef;
            Je=(Ep*Ep'+Eu*Eu'+Ef*Ef')*param.Je;
        end
        % Add hard constraints
        c=[];
        if (param.hard.Jt>0)
            c=[c;Jt-param.hard.Jt];
        end
        if (param.hard.Je>0)
            c=[c;Je-param.hard.Je];
        end
        ceq = [];
    end

Fu = (C*(t_hat)'- ysp)*C*(Wu');
H = ((Wu)*C')*(C*(Wu'));
H = H + M;

% Flow Constraints
vlb=(lower'-xmean2(np+1:1:nb)')./xstd2(np+1:1:nb)';
vub=(upper'-xmean2(np+1:1:nb)')./xstd2(np+1:1:nb)';
vlb=vlb-U_future;
vub=vub-U_future;
% Volume constraints
A=[-eye(length(vlb));eye(length(vlb))];
b=[-vlb; vub];
Av=1*ones(1,length(upper))*diag(xstd2(1,np+1:1:nb));
bv=vol_con-1*ones(1,length(upper))*xmean2(1,np+1:nb)'-1*ones(1,length(upper))*(diag(xstd2(1,np+1:1:nb))*U_future);
A=[A;Av];
b=[b;bv];

%Calculate confidence limits for Jt
sa=diag(var(T));
T2=sort(diag(T/sa*T'));
pos=round(.95*length(T2));
lambdat=1/(T2(pos));
%Calculate confidence limits for Je
spe=sort(diag(E*E',0));
pos=round(.95*length(spe));
lambdae=1/(spe(pos)+2);
%Jt,Je constraints
hard=struct('Jt',0,'Je',0,'fix',0);
param=[];
param.sa=sa;
param.tp=t_hat;
param.td=Wu;
param.Jt=lambdat;

%param.nu=(t_hat)*C';
%param.vuf=(Wu)*C';
Xc=[Xp' X2'];
param.ep=Xc*(eye(length(Xc))-Worg*Porg');
param.eu=eye(size(Wu*Pu'))-Wu*Pu';
param.ef=betau*Pf'*(eye(size(Wf*Pf'))-Wf*Pf');
param.Je=lambdae;
%% Optimisation
options = optimset('Algorithm','active-set','Display','off');%'active-set'(number of iterations exceeded),'interior-point'(change in x smaller than tol),'sqp'(no feasible solution found),'trust-region-reflective'(does not support Je Jt constraints)
options.MaxIter=10000;
options.MaxFunEvals=options.MaxIter;
optionsx0 = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
f = zeros(size(U_future)); % assumes x0 is the initial point
xnew = linprog(f,A,b,[],[],[],[],optionsx0); % Calculate a feasible x0

problem=createOptimProblem('fmincon','objective',@(x) myfun(x,H,Fu'),'x0',xnew,'Aineq',A,'bineq',b,'nonlcon',@(x) mycon(x,param),'options',options);
gs=GlobalSearch('Display','off');
[deltaU,fval,exitflag,bvoutput,lambda] =run(gs,problem);

%fprintf('Exitflag :%2.1f',exitflag); %Print exitflag
tk=param.tp+deltaU'*param.td;
U_future = U_future+deltaU;
Ufuture = U_future.*xstd2(1,np+1:1:np+nu)' + xmean2(1,np+1:1:np+nu)';

yest = tk*C';
yest = yest.*ystd+ymean;
end


