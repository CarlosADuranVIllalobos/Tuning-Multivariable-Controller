%% MPC_main_Indepensim.m
% Function implements MPC control- only runs 1 simulation per run
% Copyright Carlos Alberto Duran Villalobos August 2021 - University of Manchester.
% Last modified: March , 2022
% For indpensim_run.m from Stephen Goldrick September 2014

%% close all past data and load libraries
close all
clc
clear all
addpath('Tuning_methods');           %Add tuning functions folder
load('training_batches.mat')    %Load training dataset, batch_ini and batch_rep
rep_names=fieldnames(training); %get the names of the replicates

%% initialize variables
batch_opt = 70;                 % batches optimised
u_disturb=.20;                  % percentage max amplitude rbs disturbance
ysp=30;                         %ideal sp for penicillin
ffactor=1;                      %forgetting factor
num_x=6;                        %number of controlled variables
num_u=2;                        %number of manipulated variables
ustart=50;                      %Initial control action
invar=6;                        %inital conditions used for prediction
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
Recipe_Fs=[5:5:1150];
Recipe_Fs_nom = [8 8 8 15*ones(1,9) 30*ones(1,4) 75*ones(1,4) 150*ones(1,4) 30*ones(1,4) 37*ones(1,4) 43*ones(1,4) 47*ones(1,4) 51*ones(1,4) 57*ones(1,4) 61*ones(1,4) 65*ones(1,4) 72*ones(1,4) 76*ones(1,4) 80*ones(1,4) 84*ones(1,4) 90*ones(1,4) 136*ones(1,4) 130*ones(1,80) 120*ones(1,70)];
%Rec_Fs = [3 5 8 50*ones(1,227) ]; %Suboptimal
Recipe_PAA=[5:5:1150];  %set the sampling time to 1 hour
Recipe_PAA_nom = [5*ones(1,5)  0*ones(1,35)    10*ones(1,160)    0*ones(1,30)  ];%Nominal trajectory for FPaa
timeinterval=length(Recipe_Fs_nom);

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
    online_MVT_PAA=[]; %To observe changes in the MVT throughout the batch
    
    %% Identify initial set
    Seed_ref =Seed_ref+batch_ini;
    MVT_Fs = training.(char(rep_names(m))).MVT_Fs;
    MVT_PAA = training.(char(rep_names(m))).MVT_PAA;
    X = training.(char(rep_names(m))).X;
    Y = training.(char(rep_names(m))).Y;
    output=Y;       %Real output value if using forgetting factor
    pre_Y=Y; %pre optimised output
    
    %%% normalize
    [X2,xmean,xstd]=zscore(X);
    [Y2,ymean,ystd]=zscore(Y);
    xstd(xstd<=1e-12)=1e-12;
    
    %%% Choosing number of lv in MC-CV
    MCCV=plsmccv(X2,Y2,size(X2,1),10000);
    lv=MCCV.optLV;

    %%%Get PLS model
    [BETA,W,T,P,Q,Wp,E,F]=pls(X2,Y2,lv);
    E=X2-T*P';
    
    %%% Selecting M
    %%%% Constant approach
    %M=0*ones(length(timeinterval*u_nom),1);
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
        bt=length(Recipe_Fs_nom); %lenght of trajectories per batch
        x0=[15+0.5*randn 297+0.5*randn 6.5+0.1*randn 5.800e+04+500*randn 0.038+0.001*randn 0.20+0.05*randn]; %initial conditions
        contr_param.Recipe_Fs_sp=Recipe_Fs_nom;
        contr_param.Recipe_PAA_sp=Recipe_PAA_nom;
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
        addpath('IndPensim');
        [xref] = INDSIM2(Recipe_Fs, contr_param.Recipe_Fs_sp, Recipe_PAA, contr_param.Recipe_PAA_sp, x0, 1, Seed_ref, Rand_ref);
        [xc,yc_nom]=convfactor(xref,ustart); %convert values for PLS
        pre_Y=[pre_Y;yc_nom];
        rmpath('IndPensim')
        
        %%% run batch simulation with MPC
        addpath('IndPensim_Fast_MPC');
        [xref,unom_Fs,unom_PAA] = INDSIM2(Recipe_Fs, contr_param.Recipe_Fs_sp, Recipe_PAA, contr_param.Recipe_PAA_sp, x0, 1, Seed_ref, Rand_ref, contr_param);
        [xc,yc]=convfactor(xref,ustart); %convert values for PLS
        xc=[x0 xc];
        errest=sqrt((yc-yest(2))^2);
        
        fprintf('Batch:%3i     Nominal Y:%6.4f     Ouput:%6.4f        Errorest:%6.4f\n',[length(output)+1, yc_nom, yc, errest]);
        rmpath('IndPensim_Fast_MPC')
        %%% Update MV_trajectories and estimation error
        MVT_Fs(:,i+batch_ini) = unom_Fs;
        MVT_PAA(:,i+batch_ini) = unom_PAA;
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
%         M=M*ones(timeinterval*num_u,1);
        
        %%% M for confidence intervals approach
        %M = tune_btsp(minM, maxM, length(MVT_Fs), num_u, num_x, X2, Y2, ustart, invar, BETA, lv);
        M = tune_btsp2(minM, maxM, timeinterval, num_u, num_x, X2, Y2, ustart, invar, BETA, lv);
    end
    %% Store results 
    Fast_MPC.(char(rep_names(m))).output = output; % Output after the control action
    Fast_MPC.(char(rep_names(m))).error_end = error_end; % Estimation Error
    Fast_MPC.(char(rep_names(m))).MVT_Fs = MVT_Fs; %Final MVT
    Fast_MPC.(char(rep_names(m))).MVT_PAA = MVT_PAA; %Final MVT

end
%% Save results
FileName=['Results/BSCL2']; %Change for another save name
save(FileName,'Fast_MPC');

figure(1)
set(gcf,'name','Output quality batch run')
set(gcf, 'color', [1 1 1])
set(gca,'color', [1 1 1],'fontsize',20)
hold on;
plot(output,'black','linewidth',1.5);
set(get(gca,'YLabel'),'String','Final Penicillin concentration g/l','fontsize',20)
set(get(gca,'XLabel'),'String','Batch number','fontsize',20)
plot(pre_Y,'--red','linewidth',1.5);
x=[0 80];
y=[30 30];
line(x,y);
legend('YMPC','Ynominal','sp');
%exit



