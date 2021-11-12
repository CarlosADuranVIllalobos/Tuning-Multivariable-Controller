%% MPC_main.m
% Function implements MPC control and shows the output at each control point
% For mammalianFBLoopMPC.m 
% Copyright: Carlos Alberto Duran Villalobos August 2021 - University of Manchester.
% Last modified: November , 2021

%% close all past data and load libraries
close all
clc
clear all

addpath('CHOcells');
load('training_batches_rbs.mat') %Load training dataset, batch_ini and batch_rep
rep_names=fieldnames(training); %get the names of the replicates
MPC=struct; %Structure containing output data
%% initialize variables
batch_opt=70;  % batches optimised
ysp=2.84e9; %ideal sp for cells
ffactor=1; %forgetting factor
var_ini = 0.20;    %Initial variability
num_x=5; %number of controlled variables
num_u=1; %number of manipulated variables
ustart=45; % Initial control action
invar=5; % inital conditions used for prediction

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
    
    %normalize
    [X2,xmean,xstd]=zscore(X);
    [Y2,ymean,ystd]=zscore(Y);
    xstd(xstd<=1e-12)=1e-12;
    
    rng(Seed_ref); %set seed
    %Choosing number of lv in MC-CV
    MCCV=plsmccv(X2,Y2,size(X2,1),10000);
    lv=MCCV.optLV;
    %lv=3;
    %Get PLS model
    [BETA,W,T,P,Q,Wp,E,F]=pls(X2,Y2,lv);
    E=X2-T*P';
    
    %% optimisation loop
    
    for i=1:1:batch_opt %number of optimisation batches
        %Set conditions for the optimisation
        Rand_ref =1*m; % to make sure the same random number generator works for Indpensim
        rng(Seed_ref); %set seed
        bt=length(Rec_Fs); %lenght of trajectories per batch
        x0const=[0.33e9; 0.3e9; 0.03e9; 20; 4; 20; 1.5; 0; 7]; 
        x0_var = var_ini-2*var_ini.*rand(length(x0const),1); % +-  variation in initial conditions
        x0 = x0const + x0_var.*x0const; % add variability to initial conditions        
        opt_Fs=Rec_Fs; %Optimal feed
        
        %run batch simulation with nominal feed
        [Xref] = mammalianFBLoopMPC(x0,timeinterval,1,Rec_Fs);
        [xc,yc]=convfactor(Xref,ustart); %convert values for PLS
        xc=[x0(4) x0(5) x0(6) x0(7) x0(9) xc];
        pre_Y=[pre_Y;yc];
        fprintf('Nominal Y:%6.4f     Batch:%3i \n',[yc,length(output)+1]);
        
        %% MPC (This is the part to change to try a different controller + the function mvtoptu1)
        for cp=ustart:12:bt %controlpoints every 12h
            %Arrange for MPC
            [Xtrain,xp,Porg,Worg,Eorg,xmean2,xstd2]=orgplsmpc(X2,P,W,E,xmean,xstd,num_x,num_u,bt,ustart,xc,cp,invar);  %organize variables for MPC
            xk=[xp opt_Fs(cp:end)]; %known variables
            np=length(xp);
            nk=length(xk);
            nu=length(opt_Fs(cp:end));
            xk = (xk-xmean2(1:nk))./xstd2(1:nk);
            
            %Feed constraints
            lower=0*ones(1,timeinterval-cp+1); %Lower limit Fs
            upper=0.04*ones(1,timeinterval-cp+1); %Upper limit Fs
            
            %Begin optimisation
            M=1e-3*eye(nu); %MVT control tunning
            vol_con =15-x0(9); %constraint on volume
            ysppls=(ysp-ymean)./ystd;%
            [ufuture,H,Fu,yest,tk]=mvtoptu1(Xtrain,ysppls,xk,np,nu,nk,Eorg,T,Porg,Worg,Q,...
                xmean2,xstd2,ymean,ystd,upper,lower,M,vol_con);
            
            %Run batch simulation with the optimised MVT
            opt_Fs(cp:end)=ufuture(1:timeinterval-cp+1)';
            Rand_ref =1*m; %to make sure the same random number generator works for Indpensim
            [Xref] = mammalianFBLoopMPC(x0,timeinterval,1,opt_Fs);
            [xc,yc]=convfactor(Xref,ustart); %convert values for PLS
            xc=[x0(4) x0(5) x0(6) x0(7) x0(9) xc];

            errest=sqrt((yc-yest)^2);
            online_MVT_Fs=[online_MVT_Fs; opt_Fs];
            fprintf('Cp:%6.4f Yest:%6.4f output:%6.4f Errorest:%6.4f \n',[cp,yest,yc,errest]);
        end
        %Update MV_trajectories and estimation error
        MVT_Fs(:,i+batch_ini) = opt_Fs;
        error_end=[error_end; errest];
        output=[output;yc];
        
        %Update the PLS model
        X=[X*ffactor;xc];
        Y=[Y*ffactor;yc];
        [X2,xmean,xstd]=zscore(X);
        [Y2,ymean,ystd]=zscore(Y);
        xstd(xstd<=1e-12)=1e-12; % Avoid larger variations with a small number of batches
        MCCV=plsmccv(X2,Y2,10,10000);
        lv=MCCV.optLV;
        [BETA,W,T,P,Q,Wp,Rx,Ry]=pls(X2,Y2,lv);
        E=X2-T*P';
        Seed_ref =Seed_ref+1;
    end
    
%% Save the results in a structure    
    MPC.(char(rep_names(m))).output = output; % Output after the control action
    MPC.(char(rep_names(m))).error_end = error_end; % Estimation Error
    MPC.(char(rep_names(m))).MVT_Fs = MVT_Fs; %Final MVT
    MPC.(char(rep_names(m))).online_MVT_Fs = online_MVT_Fs; %Example of MVT change  
end

FileName=['Results/MPC_1e-3']; %Change for another save name
save(FileName,'MPC');

figure(1)
set(gcf,'name','Output quality batch run')
set(gcf, 'color', [1 1 1])
set(gca,'color', [1 1 1],'fontsize',20)
hold on;
plot(output,'black','linewidth',1.5);
set(get(gca,'YLabel'),'String','Final viable cell density [cells/L]','fontsize',20)
set(get(gca,'XLabel'),'String','Batch number','fontsize',20)
plot(pre_Y,'--red','linewidth',1.5);
x=[0 80];
y=[2.84e9 2.84e9];
line(x,y);
legend('YMPC','Ynominal','sp');




