%%INDSIM.m Version 2 
% Objective: Wrapper file for simulation
%% Copyright
% Stephen Goldrick Aug 2017
% Univeristyu College London, University of Manchester and Perceptive Engineering

% Modified by Carlos Alberto Duran Villalobos Aug 2017
% University of Manchester 

% Please reference  "The Development of an Industrial Scale Fed-Batch
% Fermentation Simulation", Stepen Goldrick, Andrei Stefen, David Lovett,
% Gary Montague, Barry Lennox, Published in Jan, Journal of Biotechnology 2015
%
% All rights reserved. Copyright (c) Newcastle University, University of Manchester and Perceptive Engineering.
% function [Xref] = indpensim_run 
function [Xref,unom_Fs,unom_PAA] = INDSIM2(Recipe_Fs,Recipe_Fs_sp,Recipe_PAA,Recipe_water_PAA,xini,Batch_no,Seed_ref,Rand_ref,contr_param)
 set(0,'DefaultFigureWindowStyle','docked') 
% close all
% clc
% clear all

%%  Set up simulation flags
Ctrl_flags.SBC= 0;      % Sequential Batch Control flag (SBC) = 0 -  Sequential Batch Control (Fixed trajectories)
Ctrl_flags.PRBS= 0;     % 1 - Add PRBS flucations to Substrate feed rate (Fs)
                        % 0 - Fixed Substrate feed rate (Fs) for each batch
Ctrl_flags.Fixed_Batch_length = 0; % 0 - Fixed Batch length
                                   % 1 - Uneven batch length
Ctrl_flags.IC = 0;      % Initial Conditions flag (IC) = 0 -  Randomly calculated initial conditions must control using SBC
Ctrl_flags.Inhib = 0;   % Inhibition  flag (Inhib)     = 0 - No inhibition
                        % Inhibition  flag (Inhib)     = 1 - Inhibition DO2, T, pH
                        % Inhibition  flag (Inhib)     = 2 - Inhibition of DO2,T,pH,CO_2_L,PAA and N
Ctrl_flags.Dis = 1;     % Disturbance flag  (Dis)      = 0 - No process disturbances
                        % Disturbance flag  (Dis)      = 1 - In batch fluctuations on mu_P, mu_x, c_s, c_oil,abc,PAA_c,Tcin and O2_in            
Ctrl_flags.Faults= 0;   % Fault flag  (Faults)         = 0 - No Faults
                        % Fault flag  (Faults)         = 1 - Aeration rate fault
                        % Fault flag  (Faults)         = 2 - Vessel back pressure  fault
                        % Fault flag  (Faults)         = 3 - Substrate feed rate fault
                        % Fault flag  (Faults)         = 4 - Base flowrate fault
                        % Fault flag  (Faults)         = 5 - Coolant flowrate fault
                        % Fault flag  (Faults)         = 6 - All of the above faults
                        % Fault flag  (Faults)         = 7 - Temperature sensor error
                        % Fault flag  (Faults)         = 8 - pH sensor error
Ctrl_flags.Vis= 0;      % Viscosity flag (Vis)          = 0 - uses simulated viscosity
                                               
% Off-line measurement sampling rate and analysis delay
Ctrl_flags.Off_line_m =  12;     % Off-line measurement sampling rate (hours)
Ctrl_flags.Off_line_delay =  4;  % Off-line measurement analysis time delay (hours)

%% Standard batch simulation with randomised initial conditions and batch
if(Ctrl_flags.IC==0)
    Ctrl_flags.SBC= 0; % Standard batch simulation must be controlled through sequential batch control
    Ctrl_flags.Vis= 0; % Standard batch simulation must use simulated viscosity
    color_number =1;
     Optimum_Batch_lenght = 230; % Standard batch length (h)
    if Ctrl_flags.Fixed_Batch_length ==1   
    Batch_length_variation = 25*randn(1);  % Batch length deviation (h)
    T = Optimum_Batch_lenght+Batch_length_variation; % Batch simulation time  (h)
    T= round(T);
    elseif Ctrl_flags.Fixed_Batch_length ==0
        T = Optimum_Batch_lenght;
    end 
    %Enbaling seed for repeatable random numbers for different batches
    %Seed_ref = 31; 
    % Rand_ref =1;
    % Defining consistent seed for random number generator for each
    % variable
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    intial_conds = 0.5+0.05*randn; % X_0 (g h^-1)
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    x0.mux = 0.41+0.025*randn;% Maximum specific growth rate of biomass [h^-1]
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    x0.mup = 0.041+ 0.0025*randn;  % maximum specific growth rate of penicillin [h^-1]
    h = 0.2; % Simulation sampling rate is 12 mins (h)
    %% Initialising simulation 
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    x0.S = 1 + 0.1*randn;      % substrate concentration [g L^{-1}]
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    %x0.DO2 = 15 + 0.5*randn  ; % initial dissolved oxygen concentration [mg L^{-1}]
    x0.DO2 = xini(1); 
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    x0.X = intial_conds + 0.1*randn; % initial biomass concentration [g L^{-1}]
    x0.P = 0;       % penicillin concentration [g L^{-1}]
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    %x0.V =  5.800e+04 + 500*randn;   % initial  volume [L]
    x0.V = xini(4); 
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    x0.Wt = 6.2e+04 + 500*randn; % Initial  Weight [Kg]
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    %x0.CO2outgas = 0.038 + 0.001*randn;   % initial carbon dioxide concentration offgas  [% of offgas]
    x0.CO2outgas = xini(5); 
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    %x0.O2 = 0.20+ 0.05*randn;% initial oxygen concentration offgas [% of offgas]
    x0.O2 = xini(6); 
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    %x0.pH = 6.5 + 0.1*randn;    % initial pH
    x0.pH = xini(3);
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    %x0.T =297+ 0.5*randn;     % temperature [K]
    x0.T = xini(2); 
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    x0.a0 = intial_conds*(1/3);   % type a0 biomass concentration [g L^{-1}]
    x0.a1 = intial_conds*(2/3);   % type a1 biomass concentration [g L^{-1}]
    x0.a3 =  0    ;   % type a3 biomass concentration [g L^{-1}]
    x0.a4 = 0 ;   % type a4 biomass concentration [g L^{-1}]
    x0.Culture_age =0; % initial culture age time A_t [hr]
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    x0.PAA = 1400+50*randn;% initial PAA conc [mg L^{-1}]
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    x0.NH3 = 1700+50*randn; % initial Nitrogen concentration [mg L^{-1}]
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    alpha_kla = 85 + 10*randn; % kla_constant [-]
    rng(Seed_ref+Batch_no +Rand_ref); 
    Rand_ref =Rand_ref+1;
    PAA_c =530000  + 20000*randn ; % Concentration of PAA in F_{PAA} [mg/L]
    rng(Seed_ref+Batch_no +Rand_ref); 
    N_conc_paa = 2*75000 + 2000*randn;%  Concentration of N in F_{PAA} [mg/L]
    Batch_time = 0:h:T;
    %% Temperature and pH Set points for Batch
    Ctrl_flags.T_sp = 298; % Temperature Set-point (K)
    Ctrl_flags.pH_sp = 6.5; % pH Set-point (-)
end
 rng(10);        
%% Creates process disturbances on growth rates as well as process inputs
% using a low pass filter
b1 = 1 - 0.995;
a1 = [1 -0.995];
% Penicillin specific growth rate disturbance: with SD of +/- 0.0009 [hr^{-1}]
v = randn(T/h+1, 1);
distMuP = filter(b1,a1,0.03*v);
Xinterp.distMuP = createChannel('Penicillin specific growth rate disturbance','g/Lh','h',Batch_time,distMuP);
% Biomass specific growth rate disturbance: with SD  +/- 0.011 [hr^{-1}]
v = randn(T/h+1, 1);
distMuX = filter(b1,a1,0.25*v);
Xinterp.distMuX = createChannel('Biomass specific  growth rate disturbance','hr^{-1}','h',Batch_time,distMuX);
% Substrate inlet concentration disturbance: +/- 15 [g L^{-1}]
v = randn(T/h+1, 1);
distcs = filter(b1,a1,300*v);
Xinterp.distcs = createChannel('Substrate concentration disturbance ',' g L^{-1}','h',Batch_time,distcs);
% Oil inlet concentration disturbance: +/- 15 [g L^{-1}]
v = randn(T/h+1, 1);
distcoil = filter(b1,a1,300*v);
Xinterp.distcoil = createChannel('Substrate concentration disturbance ',' g L^{-1}','h',Batch_time,distcoil);
% Acid/Base molar inlet concentration disturbance: +/- 0.004 [mol L^{-1}]
v = randn(T/h+1, 1);
distabc = filter(b1,a1,0.2*v);
Xinterp.distabc = createChannel('Acid/Base concentration disturbance ',' g L^{-1}','h',Batch_time,distabc);
v = randn(T/h+1, 1);
distPAA = filter(b1,a1,300000*v);
Xinterp.distPAA = createChannel('Phenylacetic acid concentration disturbance ',' g L^{-1}','h',Batch_time,distPAA);
% PAA inlet concentration disturbance: +/- 120  [g  L^{-1}]
v = randn(T/h+1, 1);
distPAA = filter(b1,a1,300000*v);
Xinterp.distPAA = createChannel('Phenylacetic acid concentration disturbance ',' g L^{-1}','h',Batch_time,distPAA);
% Coolant temperature inlet concentration disturbance: +/- 3 [K]
v = randn(T/h+1, 1);
distTcin = filter(b1,a1,100*v);
Xinterp.distTcin = createChannel('Coolant inlet temperature disturbance ','K','h',Batch_time,distTcin);
% Oxygen inlet concentration: +/- 0.009 [%]
v = randn(T/h+1, 1);
distO_2in = filter(b1,a1,0.02*v);
Xinterp.distO_2in = createChannel('Oxygen inlet concentration','%','h',Batch_time,distO_2in);

%% Import parameter list
par = Parameter_list(x0,alpha_kla,N_conc_paa,PAA_c);
%% Runs simulation
%disp('Running IndPenSim...');

[Xref,unom_Fs,unom_PAA] = indpensim(@fctrl_indpensim, Xinterp, x0, h, T,2,par,Ctrl_flags,Recipe_Fs,Recipe_Fs_sp,Recipe_PAA,Recipe_water_PAA,contr_param);





