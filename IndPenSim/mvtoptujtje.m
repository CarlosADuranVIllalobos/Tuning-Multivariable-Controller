%Function to optimise an MVc_param.T
%c_param.Qopyright c_param.Qarlos Alberto Duran Villalobos July 2017-University of Manchester

function [Ufuture,H,Fu,tk]=mvtoptujtje(Xtrain,ysp,Xc,np,nu,nb,E,c_param,Porg,Worg,xmean2,xstd2,upper,lower,M,vol_con)
global yest

% Xtrain: c_param.Training values for the PLS model
% ysp: set-point for y
% Xc: nominal vector of predictors
% np: length of past predictors
% nu: length of predictors for the MVc_param.T
% nb: length of past predictors + MVc_param.T
% E: Residuals PLS
% c_param.T: Scores PLS
% Porg: X Loadings PLS organised by control points
% Worg: Weights PLS organised by control points
% c_param.Q: Y LoadingsPLS
% xmean2: mean of predictors organised by control points
% xstd2: std of predictors organised by control points
% c_param.ymean: mean of outputs
% c_param.ystd: std of outputs
% upper: Upper constraints MVc_param.T
% lower: lower constraints MVc_param.T
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

%Mising data by c_param.Trimed Score Regression
Xpu = Xtrain(:,1:np+nu);
Xf = Xtrain(:,np+nu+1:end);
beta=Wpu/(Wpu'*Xpu'*Xpu*Wpu)*Wpu'*Xpu'*Xtrain*Worg;
t_hat=Xp'*beta;
X2=(t_hat*Pf')';
%D=beta*Pf'*Wf;
%Du=D(np+1:np+nu,:)';
%betau=beta(np+1:np+nu,:);
%Mising data by Known Data Regression
%       theta=c_param.T'*c_param.T;
%      % Beta=(Ppu*theta*Ppu')\Ppu*theta*Pf'*Wf+Wpu; %PLS
%       Beta=theta*Ppu'/(Ppu*theta*Ppu');
%       t_hat=Beta*Xp;
%Missing data by PMP
%       D =((Ppu'*Ppu)\Ppu');
%       Du=D(:,np+1:np+nu);
%     X2=Pf*D*Xp;

% Optimisation Routine tracking YSP with c_param.TSR
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
            t=param.tp+x'*param.td;
            Jt=t/param.sa*t'/param.Jt;
        end
        if (param.hard.Je>0)
            %Ep=param.ep;
            %Eu=x'*param.eu;
            %Ef=x'*param.ef;
            eu=(param.u_nom + x)'*param.eu;
            %Je=(Ep*Ep'+Eu*Eu')*param.Je;
            Je=eu*eu'*param.Je;
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

Fu = (c_param.Q*(t_hat)'- ysp)*c_param.Q*(Wu');
H = ((Wu)*c_param.Q')*(c_param.Q*(Wu'));
H = H + M;

% Flow c_param.Qonstraints
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

%Jt,Je constraints
hard=struct('Jt',c_param.jt,'Je',c_param.je);
param=[];
param.u_nom=U_future;
param.hard=hard;
param.tp=t_hat;
param.td=Wu;
c_param.Tu = Xtrain(:,np+1:np+nu)*Wu;

%calculate confidence limits for Jt
if (param.hard.Jt>0)   
    T=c_param.T(c_param.Y_ind,:);
    sa=diag(var(T));
    ev=sort(diag(T/sa*T'));
    cIt =max(bootci(1000,{@median,ev},'alpha',.01));
    lambdat=1/cIt;   
    param.Jt=lambdat;
end

%Calculate confidence limits for Je
if (param.hard.Je>0)
    U_nompast = Xtrain(c_param.Y_ind,np+1:np+nu);
    Eu=(U_nompast-c_param.Tu(c_param.Y_ind,:)*Pu');
    ev=sort(diag(Eu*Eu'));
    c_param.QIe = max(bootci(1000,{@median,ev},'alpha',.01));
    lambdae=1/c_param.QIe;
    Xc=[Xp' X2'];
    %param.ep=Xc*(eye(length(Xc))-Worg*Porg');
    param.eu=eye(size(Wu*Pu'))-Wu*Pu';
    %param.ef=betau*Pf'*(eye(size(Wf*Pf'))-Wf*Pf');
    param.Je=lambdae;
end

%% Optimisation
options = optimset('Algorithm','active-set','Display','off');%'active-set'(number of iterations exceeded),'interior-point'(change in x smaller than tol),'sqp'(no feasible solution found),'trust-region-reflective'(does not support Je Jt constraints)
options.MaxIter=10000;
options.MaxFunEvals=options.MaxIter;
optionsx0 = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
f = zeros(size(U_future)); % assumes x0 is the initial point
if ~isreal(b)
    print(b)
end

xnew = linprog(f,A,b,[],[],[],[],optionsx0); % Calculate a feasible x0

problem=createOptimProblem('fmincon','objective',@(x) myfun(x,H,Fu'),'x0',xnew,'Aineq',A,'bineq',b,'nonlcon',@(x) mycon(x,param),'options',options);
gs=GlobalSearch('Display','off');
[deltaU,fval,exitflag,bvoutput,lambda] =run(gs,problem);

%fprintf('Exitflag :%2.1f',exitflag); %Print exitflag
tk=param.tp+deltaU'*param.td;
U_future = U_future+deltaU;
Ufuture = U_future.*xstd2(1,np+1:1:np+nu)' + xmean2(1,np+1:1:np+nu)';

yest_pre = param.tp*c_param.Q';
yest_pre = yest_pre.*c_param.ystd+c_param.ymean;
% if nu==101% First cp, change depending on the simulation
if nu==362 % First cp, change depending on the simulation
    yest(1)=yest_pre;
end

yest_aft = tk*c_param.Q';
yest_aft = yest_aft.*c_param.ystd+c_param.ymean;
yest(2) = yest_aft;
end


