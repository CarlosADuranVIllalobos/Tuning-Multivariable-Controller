%% plot_results.m
% Function plot results of MPC
% For indpensim_run.m from Stephen Goldrick September 2014
% Copyright Carlos Alberto Duran Villalobos August 2017 - University of Manchester.
% Last modified: March , 2021

% IndPenSim and sub-routines are Copyright (c) Newcastle University, University of Manchester and Perceptive Engineering.
% Please reference  "The Development of an Industrial Scale Fed-Batch
% Fermentation Simulation", Stepen Goldrick, Andrei Stefen, David Lovett, Gary Montague, Barry Lennox, Submitted to the Journal of Biotechnology % September  2014

%% close all past data and load libraries
close all
clc
clear all
addpath('Results')

load('nominal_batches.mat');
rep_names = fieldnames(nominal); %get the names of the replicates
y_nom = nominal.(char(rep_names(1))).Y;

load('MPC_1e-3.mat');
y_mpc1 = MPC.(char(rep_names(1))).output;
load('FastMPC_1e-3.mat');
y_mpc2 = Fast_MPC.(char(rep_names(1))).output;


figure(1)
plot(y_mpc1,'black','linewidth',1.5);
hold on;
plot(y_mpc2,'blue','linewidth',1.5);
% plot(y_mpc4,'green','linewidth',1.5);
% plot(y_mpc5,'magenta','linewidth',1.5);
set(get(gca,'YLabel'),'String','Final viable cell density [cells/L]','fontsize',20)
set(get(gca,'XLabel'),'String','Batch number','fontsize',20)
plot(y_nom,'--red','linewidth',1.5);
x=[0 80];
y=[2.84e9 2.84e9];
line(x,y);
legend('MPC M=1e-3','FastMPC M=1e-3','Ynominal','Ysp');
set(findall(gcf,'-property','FontSize'),'FontSize',34)


figure(2)
boxes=[y_nom, y_mpc1 , y_mpc2] ;
boxplot(boxes);
hold on; 
ylabel('Final viable cell density [cells/L]')
xticklabels({'Nominal','MPC M=1e-3','FastMPC M=1e-3'})
%legend('ODE generated','ODE + measure noise','Soft sensor')
set(findall(gcf,'-property','FontSize'),'FontSize',34)
