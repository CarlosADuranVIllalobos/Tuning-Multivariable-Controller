%% create_nominal.m
% Script to create batches with nominal trajectories
% For indpensim_run.m from Stephen Goldrick September 2014
% Copyright Carlos Alberto Duran Villalobos August 2020 - University of Manchester.

%% close all past data
close all
clc
clear all
addpath('IndPensim');

%% initialize variables
end_point=[];
batch_total=80;  % initial number of batches
batch_exp=3;  % experimental batches
u_disturb=0; % percentage max amplitude rbs disturbance
var_ini=.10;  % initial variability 10%
ustart=50; % Initial control action

%% recipe settings
Recipe_Fs=[5:5:1150];
Rec_Fs = [8 8 8 15*ones(1,9) 30*ones(1,4) 75*ones(1,4) 150*ones(1,4) 30*ones(1,4) 37*ones(1,4) 43*ones(1,4) 47*ones(1,4) 51*ones(1,4) 57*ones(1,4) 61*ones(1,4) 65*ones(1,4) 72*ones(1,4) 76*ones(1,4) 80*ones(1,4) 84*ones(1,4) 90*ones(1,4) 136*ones(1,4) 130*ones(1,80) 120*ones(1,70)];
%Rec_Fs = [3 5 8 50*ones(1,227) ]; %Suboptimal
Recipe_PAA=[5:5:1150];  %set the sampling time to 1 hour
Rec_PAA = [5*ones(1,5)  0*ones(1,35)    10*ones(1,160)    0*ones(1,30)  ];%Nominal trajectory for FPaa
timeinterval=length(Rec_Fs);

%% Enbaling seed for repeatable random numbers for different batches
nominal= struct;
namelist = 'Replicate_'+string([1:batch_exp]);
for m=1:1:batch_exp
    Seed_ref=31*m;
    Rand_ref =1*m;
    X=[];
    Y=[];
    MVT_Fs=[];
    MVT_PAA=[];
    
    %% run initial set
    for i=1:1:batch_total
        rng(Seed_ref);
        %x0=[15+0.5*randn 297+0.5*randn 6.5+0.1*randn 5.800e+04+500*randn 0.038+0.001*randn 0.20+0.05*randn];
        x0=[15+15*var_ini*randn 297+297*var_ini*randn 6.5+6.5*var_ini*randn 5.800e+04+5.800e+04*var_ini*randn 0.038+0.038*var_ini*randn 0.20+0.20*var_ini*randn];
        
        Recipe_Fs_sp=Rec_Fs;
        Recipe_PAA_sp=Rec_PAA;
        
        [xref] = INDSIM2(Recipe_Fs,Recipe_Fs_sp,Recipe_PAA, Recipe_PAA_sp,x0,1,Seed_ref,Rand_ref); %run batch simulation
        [xc,yc]=convfactor(xref,ustart); %convert values for PCA and PLS
        fprintf('Plant output:%6.4f     Batch:%3i \n',[yc,length(Y)+1]);
        xc=[x0 xc];
        X=[X;xc];
        Y=[Y;yc];
        MVT_Fs(:,i) = Recipe_Fs_sp;
        MVT_PAA(:,i) = Recipe_PAA_sp;
        Seed_ref =Seed_ref+1;
    end
nominal.(namelist{m}).X=X;
nominal.(namelist{m}).Y=Y;
nominal.(namelist{m}).MVT_Fs=MVT_Fs;
nominal.(namelist{m}).MVT_PAA=MVT_PAA;   
   
end

FileName=['Results/nominal_batches'];
save(FileName,'nominal');