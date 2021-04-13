%% MPC_main.m
% Function implements MPC control and shows the output at each control point
% For indpensim_run.m from Stephen Goldrick September 2014
% Copyright: Carlos Alberto Duran Villalobos August 2017 - University of Manchester.
% Last modified: March , 2021

%% close all past data and load libraries
close all
clc
clear all

addpath('IndPensim');
load('training_batches.mat') %Load training dataset, batch_ini and batch_rep
rep_names=fieldnames(training); %get the names of the replicates
MPC=struct; %Structure containing output data
%% initialize variables
batch_opt=70;  % batches optimised
ysp=30; %ideal sp for penicillin
ffactor=1; %forgetting factor
num_x=6; %number of controlled variables
num_u=2; %number of manipulated variables
ustart=50; % Initial control action
invar=6; % inital conditions used for prediction

%% recipe settings
Recipe_Fs=[5:5:1150];
Recipe_Fs_sp = [8 8 8 15*ones(1,9) 30*ones(1,4) 75*ones(1,4) 150*ones(1,4) 30*ones(1,4) 37*ones(1,4) 43*ones(1,4) 47*ones(1,4) 51*ones(1,4) 57*ones(1,4) 61*ones(1,4) 65*ones(1,4) 72*ones(1,4) 76*ones(1,4) 80*ones(1,4) 84*ones(1,4) 90*ones(1,4) 136*ones(1,4) 130*ones(1,80) 120*ones(1,70)];
%Rec_Fs = [3 5 8 50*ones(1,227) ]; %Suboptimal
Recipe_PAA=[5:5:1150];  %set the sampling time to 1 hour
Recipe_PAA_sp = [5*ones(1,5)  0*ones(1,35)    10*ones(1,160)    0*ones(1,30)  ];%Nominal trajectory for FPaa
timeinterval=length(Recipe_Fs_sp);

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
    
    %normalize
    [X2,xmean,xstd]=zscore(X);
    [Y2,ymean,ystd]=zscore(Y);
    xstd(xstd<=1e-12)=1e-12;
    
    %Choosing number of lv in MC-CV
    MCCV=plsmccv(X2,Y2,size(X2,1));
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
        bt=length(Recipe_Fs_sp); %lenght of trajectories per batch
        x0=[15+0.5*randn 297+0.5*randn 6.5+0.1*randn 5.800e+04+500*randn 0.038+0.001*randn 0.20+0.05*randn]; %initial conditions
        unom_Fs=Recipe_Fs_sp;
        unom_PAA=Recipe_PAA_sp;
        
        %run batch simulation with nominal feed
        [xref] = INDSIM2(Recipe_Fs,unom_Fs,Recipe_PAA, unom_PAA,x0,1,Seed_ref,Rand_ref);
        [xc,yc]=convfactor(xref,ustart); %convert values for PLS
        xc=[x0 xc];
        pre_Y=[pre_Y;yc];
        fprintf('Nominal Y:%6.4f     Batch:%3i \n',[yc,length(output)+1]);
        
        %% MPC (This is the part to change to try a different controller + the function mvtoptu1)
        for cp=50:24:230 %controlpoints every 24h
            %Arrange for MPC
            [Xtrain,xp,Porg,Worg,Eorg,xmean2,xstd2]=orgplsmpc(X2,P,W,E,xmean,xstd,num_x,num_u,bt,ustart,xc,cp,invar);  %organize variables for MPC
            xk=[xp unom_Fs(cp:end) unom_PAA(cp:end) ]; %known variables
            np=length(xp);
            nk=length(xk);
            nu=length(unom_Fs(cp:end))+length(unom_PAA(cp:end));
            xk = (xk-xmean2(1:nk))./xstd2(1:nk);
            
            %Feed constraints
            lower_Fs=5*ones(1,timeinterval-cp+1); %Lower limit Fs
            upper_Fs=180*ones(1,timeinterval-cp+1); %Upper limit Fs
            lower_PAA=0.01*ones(1,timeinterval-cp+1); %Lower limit PAA
            upper_PAA=15*ones(1,timeinterval-cp+1); %Upper limit PAA
            lower=[lower_Fs lower_PAA];
            upper=[upper_Fs upper_PAA];
            
            %Begin optimisation
            M=1e-3*eye(nu); %MVT control tunning
            vol_con=29000; %constraint on volume
            ysppls=(ysp-ymean)./ystd;%
            [ufuture,H,Fu,yest,tk]=mvtoptu1(Xtrain,ysppls,xk,np,nu,nk,Eorg,T,Porg,Worg,Q,...
                xmean2,xstd2,ymean,ystd,upper,lower,M,vol_con);
            
            %Run batch simulation with the optimised MVT
            unom_Fs(cp:end)=ufuture(1:timeinterval-cp+1)';
            unom_PAA(cp:end)= ufuture(timeinterval-cp+2:end)';
            Rand_ref =1*m; %to make sure the same random number generator works for Indpensim
            [xref] = INDSIM2(Recipe_Fs,unom_Fs,Recipe_PAA,unom_PAA,x0,1,Seed_ref,Rand_ref); %run batch simulation
            [xc,yc]=convfactor(xref,ustart); %convert values for PCA and PLS
            xc=[x0 xc];
            errest=sqrt((yc-yest)^2);
            online_MVT_Fs=[online_MVT_Fs; unom_Fs];
            online_MVT_PAA=[online_MVT_PAA;  unom_PAA];
            fprintf('Cp:%6.4f Yest:%6.4f output:%6.4f Errorest:%6.4f \n',[cp,yest,yc,errest]);
        end
        %Update MV_trajectories and estimation error
        MVT_Fs(:,i+batch_ini) = unom_Fs;
        MVT_PAA(:,i+batch_ini) = unom_PAA;
        error_end=[error_end; errest];
        output=[output;yc];
        
        %Update the PLS model
        X=[X*ffactor;xc];
        Y=[Y*ffactor;yc];
        [X2,xmean,xstd]=zscore(X);
        [Y2,ymean,ystd]=zscore(Y);
        xstd(xstd<=1e-12)=1e-12; % Avoid larger variations with a small number of batches
        MCCV=plsmccv(X2,Y2,10);
        lv=MCCV.optLV;
        [BETA,W,T,P,Q,Wp,Rx,Ry]=pls(X2,Y2,lv);
        E=X2-T*P';
        Seed_ref =Seed_ref+1;
    end
    
%% Save the results in a structure    
    MPC.(char(rep_names(m))).output = output; % Output after the control action
    MPC.(char(rep_names(m))).error_end = error_end; % Estimation Error
    MPC.(char(rep_names(m))).MVT_Fs = MVT_Fs; %Final MVT
    MPC.(char(rep_names(m))).MVT_PAA = MVT_PAA; %Final MVT
    MPC.(char(rep_names(m))).online_MVT_Fs = online_MVT_Fs; %Example of MVT change
    MPC.(char(rep_names(m))).online_MVT_PAA = online_MVT_PAA; %Example of MVT change  
end

FileName=['Results/MPC_1e-3']; %Change for another save name
save(FileName,'MPC');

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




