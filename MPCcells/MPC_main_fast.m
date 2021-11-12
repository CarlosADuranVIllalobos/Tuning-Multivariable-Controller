%% MPC_main_fast.m
% Function implements MPC control- only runs 1 simulation per run
% For mammalianFBLoopMPC.m 
% Copyright: Carlos Alberto Duran Villalobos August 2021 - University of Manchester.
% Last modified: November , 2021

%% close all past data and load libraries
close all
clc
clear all

load('training_batches_rbs.mat')    %Load training dataset, batch_ini and batch_rep
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
contr_param=struct;             %Structure Parameters for the controller
Fast_MPC=struct;                %Structure containing output data
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

    %Get PLS model
    [BETA,W,T,P,Q,Wp,E,F]=pls(X2,Y2,lv);
    E=X2-T*P';
    
    %No approach
    M = 1e-3*ones(length(BETA),1);
     %M for confidence intervals approach
%       n=10000;
%      mina=0;
%      maxa=1;
%      lenx=size(X2,2);
%     [x,bootsamp] = bootstrp(n,@(bootr)pls(bootr(:,1:lenx),bootr(:,lenx+1),lv),[X2,Y2]);
%     se=0;
%     for i=1:1:n
%         se=se+(BETA'-x(i,:)).^2;
%     end
%     see1=sqrt(se/(n-1));
%     see=see1*2;  %2 times std.
%     Bbt=abs(BETA)'-see;
%     indmin=find(Bbt>=0);
%     ind=find(not(Bbt>=0));
%     M=ones(length(BETA),1);
%     M(ind)= normalize(abs(Bbt(ind)), 'range', [mina maxa]);
%     M(indmin)= mina;
%     
%     %Good confidance intervals only 
%     n=10000;
%     mina=1e-12;
%     maxa=1;
%     lenx=size(X2,2);
%     [x,bootsamp] = bootstrp(n,@(bootr)pls(bootr(:,1:lenx),bootr(:,lenx+1),lv),[X2,Y2]);
%     se=0;
%     for i=1:1:n
%         se=se+(BETA'-x(i,:)).^2;
%     end
%     see1=sqrt(se/(n-1));
%     see=see1*3;  %3 times std.
%     Bbt=abs(BETA)'-see;
%     indmin=find(Bbt>=0);
%     ind=find(not(Bbt>=0));
%     M=ones(length(BETA),1);
%     M(ind)= 1e6;
%     M(indmin)= normalize(abs(Bbt(indmin)), 'range', [mina maxa]);

     %Less constrained at the beggining and more towards the end
%     lenx=size(X2,2); 
     %Exponential
%       n=10000;
%      mina=0;
%      maxa=1;
%      lenx=size(X2,2);
%      n = length(MVT_Fs)-ustart+1;
%      mina = 1e-6;
%      maxa = 1;
%      gf = (maxa/mina)^(1/(n - 1));
%      k = 1:n;
%      y = mina*gf.^(k-1);
%      M=ones(lenx,1);
%      for i=1:num_u
%          val=1;
%          for k=ustart:length(MVT_Fs)
%              index=(k-1)*num_x + (k-ustart)*num_u + i + num_x+invar;
%              M(index)=y(val);
%              val=val+1;
%          end
%      end
     
     
     %Linear
%      n = length(MVT_Fs)-ustart+1;
%      mina = 1e-12;
%      maxa = 1;
%      gf = (maxa-mina)/n;
%      k = 1:n;
%      y = mina+gf*(k-1);
%      M=ones(lenx,1);
%      for i=1:num_u
%          val=1;
%          for k=ustart:length(MVT_Fs)
%              index=(k-1)*num_x + (k-ustart)*num_u + i + num_x+invar;
%              M(index)=y(val);
%              val=val+1;
%          end
%      end
     
     
     
    %% optimisation loop
    
    for i=1:1:batch_opt %number of optimisation batches
        %Set conditions for the optimisation
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
        
        %run batch simulation with nominal feed
        addpath('CHOcells');
        [Xref] = mammalianFBLoopMPC(x0,timeinterval,1,Rec_Fs);
        [xc,yc_nom]=convfactor(Xref,ustart); %convert values for PLS
        xc=[x0(4) x0(5) x0(6) x0(7) x0(9) xc];
        pre_Y=[pre_Y;yc_nom];
        rmpath('CHOcells')
        
        %run batch simulation with MPC
        addpath('CHOcells_Fast');
        [xref,unom_Fs,yest] = mammalianFBLoopMPC(x0,timeinterval,1,Rec_Fs,contr_param);
        [xc,yc]=convfactor(xref,ustart); %convert values for PLS
        xc=[x0(4) x0(5) x0(6) x0(7) x0(9) xc];
        errest=sqrt((yc-yest)^2);
        fprintf('Batch:%3i     Nominal Y:%6.4f     Ouput:%6.4f        Errorest:%6.4f\n',[length(output)+1, yc_nom, yc, errest]);
        rmpath('CHOcells_Fast')
        %Update MV_trajectories and estimation error
        MVT_Fs(:,i+batch_ini) = unom_Fs;
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
        
        %M for confidence intervals approach  
%          [x,bootsamp] = bootstrp(n,@(bootr)pls(bootr(:,1:lenx),bootr(:,lenx+1),lv),[X2,Y2]);
%          se=0;
%          for u=1:1:n
%              se=se+(BETA'-x(u,:)).^2;
%          end
%          see1=sqrt(se/(n-1));
%          see=see1*2;
%          Bbt=abs(BETA)'-see;
%          indmin=find(Bbt>=0);
%          ind=find(not(Bbt>=0));
%          M=ones(length(BETA),1);
%          M(ind)= normalize(abs(Bbt(ind)), 'range', [mina maxa]);
%          M(indmin)= mina;          
         
%          %Second approach
%          [x,bootsamp] = bootstrp(n,@(bootr)pls(bootr(:,1:lenx),bootr(:,lenx+1),lv),[X2,Y2]);
%          se=0;
%          for u=1:1:n
%              se=se+(BETA'-x(u,:)).^2;
%          end
%          see1=sqrt(se/(n-1));
%          see=see1*3;  %3 times std.
%          Bbt=abs(BETA)'-see;
%          indmin=find(Bbt>=0);
%          ind=find(not(Bbt>=0));
%          M=ones(length(BETA),1);
%          M(ind)= 1e6;
%          M(indmin)= normalize(abs(Bbt(indmin)), 'range', [mina maxa]);
         
         Seed_ref =Seed_ref+1; 
    end
    Fast_MPC.(char(rep_names(m))).output = output; % Output after the control action
    Fast_MPC.(char(rep_names(m))).error_end = error_end; % Estimation Error
    Fast_MPC.(char(rep_names(m))).MVT_Fs = MVT_Fs; %Final MVT

end

FileName=['Results/FastMPC_1e-3']; %Change for another save name
save(FileName,'Fast_MPC');


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




