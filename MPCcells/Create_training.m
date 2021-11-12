%% Create_training.m
% Script create training batches
% For indpensim_run.m from Stephen Goldrick September 2014
% Copyright Carlos Alberto Duran Villalobos August 2017 - University of Manchester.

% IndPenSim and sub-routines are Copyright (c) Newcastle University, University of Manchester and Perceptive Engineering.
% Please reference  "The Development of an Industrial Scale Fed-Batch
% Fermentation Simulation", Stepen Goldrick, Andrei Stefen, David Lovett, Gary Montague, Barry Lennox, Submitted to the Journal of Biotechnology % September  2014

%% close all past data
close all
clc
clear all
addpath('CHOcells');

%% initialize variables
end_point=[];
batch_ini=10;  % initial number of batches
batch_rep=1;  % Number of replicates
var_ini = 0.20;    %Initial variability
u_disturb=0.02; % percentage max amplitude rbs disturbance
cp=45; % Initial control action
ustart=cp;

%% recipe settings
Rec_Fs = [zeros(1,44) 0.02*ones(1,101)];
timeinterval=length(Rec_Fs);

%% Enbaling seed for repeatable random numbers for different batches
Seed_ref = 31;
Rand_ref =1;
training= struct;
namelist = 'Replicate_'+string([1:batch_rep]);
for m=1:1:batch_rep
    Seed_ref=31*m;
    Rand_ref =1*m;
    X=[];
    Y=[];
    MVT_Fs=[];
    MVT_PAA=[];
    %% run initial set
    for i=1:1:batch_ini
        rng(Seed_ref);
        x0const=[0.33e9; 0.3e9; 0.03e9; 20; 4; 20; 1.5; 0; 7];
        x0_var = var_ini-2*var_ini.*rand(length(x0const),1); % +-  variation in initial conditions
        x0 = x0const + x0_var.*x0const; % add variability to initial conditions
        Recipe_Fs_sp=Rec_Fs;
        
        Recipe_Fs_sp(cp:end)=Recipe_Fs_sp(cp:end)+idinput(timeinterval-cp+1,'rbs',[0 .05],[-u_disturb*max(Rec_Fs) u_disturb*max(Rec_Fs)])';
   
        Recipe_Fs_sp(Recipe_Fs_sp<0)=0; %Replace with 0 if negative

        [Xref] = mammalianFBLoopMPC(x0,timeinterval,1,Recipe_Fs_sp);
        Seed_ref =Seed_ref+1;
        
        [xc,yc]=convfactor(Xref,ustart); %convert values for PCA and PLS
        fprintf('Plant output:%6.4f     Batch:%3i \n',[yc,length(Y)+1]);
        xc=[x0(4) x0(5) x0(6) x0(7) x0(9) xc];
        X=[X;xc];
        Y=[Y;yc];
        MVT_Fs(:,i) = Recipe_Fs_sp;
    end
    training.(namelist{m}).X=X;
    training.(namelist{m}).Y=Y;
    training.(namelist{m}).MVT_Fs=MVT_Fs;
end

FileName=['training_batches_rbs'];
save(FileName,'training','batch_ini','batch_rep');