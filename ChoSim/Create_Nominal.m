%% create_nominal.m
% Script to create batches with nominal trajectories
% For indpensim_run.m from Stephen Goldrick September 2014
% Copyright Carlos Alberto Duran Villalobos August 2020 - University of Manchester.

%% close all past data
close all
clc
clear all
addpath('CHOcells');

%% initialize variables
end_point=[];
batch_total=80;  % initial number of batches
batch_exp=3;  % experimental batches
var_ini = 0.20;    %Initial variability
u_disturb=0; % percentage max amplitude rbs disturbance
ustart=45; % Initial control action

%% recipe settings
Rec_Fs = [zeros(1,44) 0.02*ones(1,101)];
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
    
    %% run initial set
    for i=1:1:batch_total
        rng(Seed_ref);
        x0const=[0.33e9; 0.3e9; 0.03e9; 20; 4; 20; 1.5; 0; 7];
        x0_var = var_ini-2*var_ini.*rand(length(x0const),1); % +-  variation in initial conditions
        x0 = x0const + x0_var.*x0const; % add variability to initial conditions
        
        Recipe_Fs_sp=Rec_Fs;
        
        [xref] = mammalianFBLoopMPC(x0,timeinterval,1,Recipe_Fs_sp);
        [xc,yc]=convfactor(xref,ustart); %convert values for PCA and PLS
        fprintf('Plant output:%6.4f     Batch:%3i \n',[yc,length(Y)+1]);
         xc=[x0(4) x0(5) x0(6) x0(7) x0(9) xc];
        X=[X;xc];
        Y=[Y;yc];
        MVT_Fs(:,i) = Recipe_Fs_sp;
        Seed_ref =Seed_ref+1;
    end
nominal.(namelist{m}).X=X;
nominal.(namelist{m}).Y=Y;
nominal.(namelist{m}).MVT_Fs=MVT_Fs; 
   
end

FileName=['Results/nominal_batches'];
save(FileName,'nominal');