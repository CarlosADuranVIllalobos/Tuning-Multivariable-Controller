%% MPC_main.m
% Function implements MPC control- only runs 1 simulation per run
% For mammalianFBLoopMPC.m 
% Copyright: Carlos Alberto Duran Villalobos August 2021 - University of Manchester.
% Last modified: March , 2022

%% close all past data and load libraries
close all
clc
clear all
addpath('Tuning_methods');           %Add tuning functions folder
load('training_batches.mat')    %Load training dataset, batch_ini and batch_rep
rep_names=fieldnames(training); %get the names of the replicates

%% initialize variables
batch_opt=70;                   % batches optimised
u_disturb=.20;                  % percentage max amplitude rbs disturbance
ysp=2.84e9;                         %ideal sp for penicillin
ffactor=1;                      %forgetting factor
var_ini = 0.20;    %Initial variability
num_x=5;                        %number of controlled variables
num_u=1;                        %number of manipulated variables
ustart=45;                      %Initial control action
invar=5;                        %inital conditions used for prediction
step=12;                        %Step of the controller [h]
jt=0;                           %Hotellin's constraint- 0:No validity constraints 1: Validity constraints on Delta_u only
je=0;                           %Q constraint- 0:No validity constraints 1: Validity constraints on Delta_u only
minM = 1e-100;                  %Minimum value of M~=0 to avoid division by 0
maxM = 1e-1;                    %Maximum value of M 
contr_param=struct;             %Structure Parameters for the controller
Fast_MPC=struct;                %Structure containing output data
global yest
yest=[0,0];
%% recipe settings
Rec_Fs = [zeros(1,44) 0.02*ones(1,101)];%Suboptimal
timeinterval=length(Rec_Fs);

%% Enbaling seed for repeatable random numbers for different batches
Seed_ref = 31;
Rand_ref =1;

for m=1:1:batch_rep %number of replicates
    fprintf('Replicate:%3i \n',m);
    Seed_ref=31*m;
    X=[];
    Y=[];
    error_end=[];   %Final estimation error
    online_MVT_Fs=[]; %To observe changes in the MVT throughout the batch
    
    %% Identify initial set
    Seed_ref =Seed_ref+batch_ini;
    MVT_Fs = training.(char(rep_names(m))).MVT_Fs;
    X = training.(char(rep_names(m))).X;
    Y = training.(char(rep_names(m))).Y;
    output=Y;       %Real output value if using forgetting factor
    pre_Y=Y; %pre optimised output
    
    %%% normalize
    [X2,xmean,xstd]=zscore(X);
    [Y2,ymean,ystd]=zscore(Y);
    xstd(xstd<=1e-12)=1e-12;
    
    rng(Seed_ref); %set seed
    %%% Choosing number of lv in MC-CV
    MCCV=plsmccv(X2,Y2,size(X2,1),10000);
    lv=MCCV.optLV;

    %%% Get PLS model
    [BETA,W,T,P,Q,Wp,E,F]=pls(X2,Y2,lv);
    E=X2-T*P';
    
    %%% Selecting M
    %%%% Constant approach
    %M=0*ones(timeinterval*num_u,1);
    %%%% Differential evolution approach
    %[M,GC,CR,F,nd,popsize,pop,bests,bests_ind,count,pop_c,flag_s,cost,new_cost,trial,new_pop] = init_de(minM, maxM, timeinterval, num_u, batch_opt);
    %%%% Linear approach
    %M = tune_linear(minM, maxM, timeinterval, num_u);   
    %%%% Exponential approach
    %M = tune_exp(minM, maxM, timeinterval, num_u);
    %%%% Confidence limits approach
    %M = tune_btsp(minM, maxM, timeinterval, num_u, num_x, X2, Y2, ustart, invar, BETA, lv);
     M = tune_btsp2(minM, maxM, timeinterval, num_u, num_x, X2, Y2, ustart, invar, BETA, lv);
    %% optimisation loop
    for i=1:1:batch_opt %number of optimisation batches
        %%% Set conditions for the optimisation
        Rand_ref =1*m; % to make sure the same random number generator works for Indpensim
        rng(Seed_ref); %set seed
        bt=length(Rec_Fs); %lenght of trajectories per batch
        x0const=[0.33e9; 0.3e9; 0.03e9; 20; 4; 20; 1.5; 0; 7]; 
        x0_var = var_ini-2*var_ini.*rand(length(x0const),1); % +-  variation in initial conditions
        x0 = x0const + x0_var.*x0const; % add variability to initial conditions    
        contr_param.Recipe_Fs_sp=Rec_Fs;
        contr_param.X2=X2;
        contr_param.Y2=Y2;
        contr_param.T=T;
        contr_param.P=P;
        contr_param.Q=Q;
        contr_param.W=W;
        contr_param.E=E;
        contr_param.xmean=xmean;
        contr_param.xstd=xstd;
        contr_param.num_x=num_x;
        contr_param.num_u=num_u;
        contr_param.ustart=ustart;
        contr_param.invar=invar;
        contr_param.step=step;
        contr_param.timeinterval=timeinterval;
        contr_param.x0=x0;
        contr_param.bt=bt;
        contr_param.ysp=ysp;
        contr_param.ystd=ystd;
        contr_param.ymean=ymean;
        contr_param.M=M;
        contr_param.jt=jt;
        contr_param.je=je;
        
        %%% run batch simulation with nominal feed
        addpath('CHOcells');
        [Xref] = mammalianFBLoopMPC(x0,timeinterval,1,Rec_Fs);
        [xc,yc_nom]=convfactor(Xref,ustart); %convert values for PLS
        xc=[x0(4) x0(5) x0(6) x0(7) x0(9) xc];
        pre_Y=[pre_Y;yc_nom];
        rmpath('CHOcells')
        
        %%% run batch simulation with MPC
        addpath('CHOcells_Fast_MPC');
        [xref,unom_Fs] = mammalianFBLoopMPC(x0,timeinterval,1,Rec_Fs,contr_param);
        [xc,yc]=convfactor(xref,ustart); %convert values for PLS
        xc=[x0(4) x0(5) x0(6) x0(7) x0(9) xc];
        errest=abs(yc-yest(2))/yc;
        fprintf('Batch:%3i     Nominal Y:%6.4e     Ouput:%6.4e        Errorest:%6.4f\n',[length(output)+1, yc_nom, yc, errest]);
        rmpath('CHOcells_Fast_MPC')
        
        %%% Update MV_trajectories and estimation error
        MVT_Fs(:,i+batch_ini) = unom_Fs;
        error_end=[error_end; errest];
        output=[output;yc];
        
        %%% Update the PLS model
        X=[X*ffactor;xc];
        Y=[Y*ffactor;yc];
        [X2,xmean,xstd]=zscore(X);
        [Y2,ymean,ystd]=zscore(Y);
        xstd(xstd<=1e-12)=1e-12; % Avoid larger variations with a small number of batches
        MCCV=plsmccv(X2,Y2,10,10000);
        lv=MCCV.optLV;
        [BETA,W,T,P,Q,Wp,Rx,Ry]=pls(X2,Y2,lv); 
        E=X2-T*P';
        Seed_ref =Seed_ref+1; %Increase the counrter for the random seed
        
        %%% M for differencial evolution approach
%         dy=(abs(yest(1)-ysp)-abs(yc-ysp))/abs(yest(1)-ysp); %Calculate fitness
%         [M, flag_s, pop_c, pop, trial, new_pop, cost, bests, bests_ind] = tune_de(dy, flag_s, pop_c, popsize, CR, F, pop, trial, new_pop, cost, bests, bests_ind, nd);
%          M=M*ones(timeinterval*num_u,1);
        
        %%% M for confidence intervals approach
        %M = tune_btsp(minM, maxM, timeinterval, num_u, num_x, X2, Y2, ustart, invar, BETA, lv);
        M = tune_btsp2(minM, maxM, timeinterval, num_u, num_x, X2, Y2, ustart, invar, BETA, lv);
        
    end
    %% Store results 
    Fast_MPC.(char(rep_names(m))).output = output; % Output after the control action
    Fast_MPC.(char(rep_names(m))).error_end = error_end; % Estimation Error
    Fast_MPC.(char(rep_names(m))).MVT_Fs = MVT_Fs; %Final MVT

end
%% Save results
FileName=['Results/BSCI2']; %Change for another save name
save(FileName,'Fast_MPC');


figure(1)
set(gcf,'name','Output quality batch run')
set(gcf, 'color', [1 1 1])
set(gca,'color', [1 1 1],'fontsize',20)
hold on;
plot(output,':oblack','linewidth',1.5);
set(get(gca,'YLabel'),'String','Final viable cell density [cells/L]','fontsize',20)
set(get(gca,'XLabel'),'String','Batch number','fontsize',20)
plot(pre_Y,':sred','linewidth',1.5);
x=[0 80];
y=[2.84e9 2.84e9];
line(x,y);
legend('YMPC','Ynominal','sp');
%exit



