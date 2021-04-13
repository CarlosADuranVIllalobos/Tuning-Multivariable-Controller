function par = Parameter_list(x0,alpha_kla,N_conc_paa,PAA_c)
%% Define model parameters
%% Penicillin model parameters
mu_p = x0.mup;          % hr^{-1}           par(1)
mux_max = x0.mux;       % hr^{-1}           par(2)
ratio_mu_e_mu_b = 0.4; % (-)               par(3)
P_std_dev = 0.0015;     % g  L^{-1}         par(4)
mean_P = 0.002;         % g  L^{-1}         par(5)
mu_v = 1.71e-4;         % g g^{-1} hr^{-1}  par(6)
mu_a = (3.5e-3);        % hr^{-1}           par(7)
mu_diff =(5.36e-3);     % hr^{-1}           par(8)
beta_1  = 0.006;        % (-)               par(9)
K_b   = 0.05 ;          % g L^{-1}          par(10)
K_diff =0.75;           % g L^{-1}          par(11)
K_diff_L = 0.09;        % g L^{-1}          par(12)
K_e = 0.009;            % g L^{-1}          par(13)
K_v = 0.05;             % g L^{-1}          par(14)
delta_r = 0.75e-004;    % cm                par(15)
k_v = 3.22e-5;          % cm hr^{-1}        par(16)
D = (2.66e-11);         % cm^{2} hr^{-1}    par(17)
rho_a0 = 0.35;          % g cm^-3           par(18)
rho_d = 0.18;           % g cm^3            par(19)
mu_h = 0.003;           % hr^{-1}           par(20)
r_0 = (1.5e-4);         % cm                par(21)
delta_0 = (1e-4);       % cm                par(22)

%% Process related parameters
Y_sx = 1.85;            % g g^{-1}          par(23)
Y_sP = 0.9;             % g g^{-1}          par(24)
m_s = 0.029;            % g g^{-1} hr^{-1}  par(25)
c_oil = 1000;           % g L^{-1}          par(26)
c_s = 600;              % g L^{-1}          par(27)
Y_O2_X = 650;           % mg g^{-1}         par(28)
Y_O2_P =160;            % mg g^{-1}         par(29)
m_O2_X =  17.5;         % mg g^{-1}         par(30)
alpha_kla;              % (-)               par(31)
a =  0.38;              % (-)               par(32)
b = 0.34;               % (-)               par(33)
c= -0.38;               % (-)               par(34)
d= 0.25;                % (-)               par(35)
Henrys_c = 0.0251;      % bar L mg^{-1}     par(36)
n_imp = 3;              % (-)               par(37)
r = 2.1;                % (m)               par(38)
r_imp = 0.85;           % (m)               par(39)
Po = 5;                 % (-)               par(40)
epsilon = 0.1;          % (-)               par(41)
g = 9.81;               % m s^{-1}          par(42)
R = 8.314;              %J K^{-1} mol^{-1}  par(43)
X_crit_DO2 = 0.1;       %  (%)              par(44)
P_crit_DO2 = 0.3;       %  (%)              par(45)
A_inhib = 1;            %   -               par(46)
Tf = 288;               %  K                par(47)
Tw = 288;               %  K                par(48)
Tcin = 285;             %   K               par(49)
Th = 333;               %  K                par(50)
Tair = 290;             %  K                par(51)
C_ps = 5.9;            % kJ kg^{-1} K^{-1}  par(52)
C_pw  = 4.18;          % kJ g^{-1} K^{-1}   par(53)
dealta_H_evap = 2430.7;% kJ kg^{-1}         par(54)
U_jacket  = 36;        % kW m^{2} K^{-1}    par(55)
A_c = 105;              % m^2               par(56)
Eg = 1.488*(10^4);       % J mol^{-1}        par(57)
Ed =  1.7325*(10^5);      % J mol^{-1}        par(58)
k_g = 450;              % (-)               par(59)
k_d = 0.25*10^30;       %  (-)              par(60)
Y_QX = 25;              % kJ/g              par(61)
abc = 0.033;           % mol L^{-1}         par(62)
gamma1=(0.0325e-5);     % mol H+/g          par(63)
gamma2 =  2.5*(1.e-11); % mol H+/g          par(64)
m_ph = 0.0025;          % mol H+/g  hr^{-1} par(65)
K1=1e-5;                % mol L^{-1}        par(66)
K2=2.5e-8;              % mol L^{-1}        par(67)
N_conc_oil = 20000;     %  mg  L^{-1}        par(68)
N_conc_paa;             %  g  L^{-1}        par(69)
N_conc_shot = 400000; %  mg kg^{-1}         par(70)
Y_NX =  10;             % mg of N_2 per gram of X DW par(71)
Y_NP = 80;              % mg of N_2 per gram of X DW  par(72)
m_N = 0.03;             % g L^{-1}  hr^{-1}  par(73)
X_crit_N = 150;         % mg L^{-1}          par(74)
PAA_c;                   % mg L^{-1}         par(75)
Y_PAA_P  =  187.5;      % mg/ g P            par(76)
Y_PAA_X =  37.5000*1.2;     % mg/ g X            par(77)
m_PAA = 1.05;          % g g^{-1} P  h^{-1}   par(78)
X_crit_PAA= 24000;           % mg/L           par(79)
P_crit_PAA = 200;           % mg/L           par(80)
B_1 = (-0.6429*(10^2)); % (-)                par(81)
B_2 = -0.1825*(10^1);   % (-)                par(82)
B_3 = 0.3649;           % (-)                par(83)
B_4 = 0.1280;           % (-)                par(84)
B_5  = -4.9496e-04;     % (-)                par(85)
delta_c_o = 0.89;     % (-)                par(86)
k_3 = 0.005;            % (-)                par(87)
k1 = 0.001;             % (-)                par(88)
k2 = 0.0001;            % (-)                par(89)
t1 = 1;                 % h                  par(90)
t2 = 250;               % h                  par(91)
q_co2 = 0.123*1.1;      % mg CO_2 (g X)^{-1} par(92)
X_crit_CO2 = 7570;        %  mg L^{-1}        par(93)
alpha_evp = (5.2400e-04);  % L hr^{-1}          par(94)
beta_T  = 2.88;         % -                  par(95)
pho_g = 1.54*1000;      % kg m^{-3}          par(96)
pho_oil = 0.90*1000;    % kg m^{-3}          par(97)
pho_w = 1000;           % kg m^{-3}          par(98)
pho_paa = 1000;         % kg m^{-3}          par(99)
O_2_in = 0.21;          % (%)                par(100)
N2_in = 0.79;           % (%)                par(101)
C_CO2_in = 0.033;       % (%)                par(102)
Tv = 373;               % K                  par(103)
T0 = 273;               % K                  par(104)
alpha_1 =  2451.8; % kJ/m^3                  par(105)



par = [mu_p,mux_max,ratio_mu_e_mu_b,P_std_dev,mean_P,mu_v,mu_a,mu_diff,... 
    beta_1,K_b,K_diff,K_diff_L,K_e,K_v,delta_r,k_v,D,rho_a0,rho_d,mu_h,r_0,delta_0... % 20 
    Y_sx,Y_sP,m_s,c_oil,c_s,Y_O2_X,Y_O2_P,m_O2_X,alpha_kla,a,b,c,d,... % 33
    Henrys_c,n_imp,r,r_imp,Po,epsilon,g,R,X_crit_DO2,P_crit_DO2,A_inhib,...
    Tf,Tw,Tcin,Th,Tair,C_ps,C_pw,dealta_H_evap, U_jacket,A_c,Eg,Ed,k_g,...
    k_d,Y_QX,abc,gamma1,gamma2,m_ph,K1,K2,N_conc_oil,N_conc_paa,...
    N_conc_shot,Y_NX,Y_NP,m_N,X_crit_N,PAA_c,...
Y_PAA_P,Y_PAA_X,m_PAA,X_crit_PAA,P_crit_PAA,B_1...
,B_2,B_3,B_4,B_5,delta_c_o,k_3,k1,k2,t1,t2,q_co2,...
X_crit_CO2,alpha_evp,beta_T,pho_g,pho_oil,pho_w,pho_paa,O_2_in,N2_in,C_CO2_in,Tv,T0,alpha_1];

%% Copyright 
% Stephen Goldrick  (c),  September, 2014
% Newcastle University/Manchester University/ Perceptive Engineering 
%
% All rights reserved. Copyright (c) Newcastle University, Manchester University and Perceptive Engineering.
%


