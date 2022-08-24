%% plot_results.m
% Function plot results of MPC
% For indpensim_run.m from Stephen Goldrick September 2014
% Copyright Carlos Alberto Duran Villalobos August 2017 - University of Manchester.
% Last modified: March , 2022

% IndPenSim and sub-routines are Copyright (c) Newcastle University, University of Manchester and Perceptive Engineering.
% Please reference  "The Development of an Industrial Scale Fed-Batch
% Fermentation Simulation", Stepen Goldrick, Andrei Stefen, David Lovett, Gary Montague, Barry Lennox, Submitted to the Journal of Biotechnology % September  2014

%% close all past data and load libraries
close all
clc
clear all

ysp=2.84e9;    


%Color palette
b=[0, 0.4470, 0.7410];
o=[0.8500, 0.3250, 0.0980];
y=[0.9290, 0.6940, 0.1250];
p=[0.4940, 0.1840, 0.5560];
g=[0.4660, 0.6740, 0.1880];
lightb= 	[0.3010, 0.7450, 0.9330];
r=[0.6350, 0.0780, 0.1840];
b2= [0, 0, 1];
g2=[0, 0.5, 0];

load('nominal_batches.mat');
rep_names = fieldnames(nominal); %get the names of the replicates
y_nom = [];
for i=1:1:length(rep_names)
    y_nom = [y_nom; nominal.(char(rep_names(i))).Y];  
end

load('1.mat');
y_mpc1 = [];
for i=1:1:length(rep_names)
    y_mpc1 = [y_mpc1; Fast_MPC.(char(rep_names(i))).output];  
end

load('e_1.mat');
y_mpc2 = [];
for i=1:1:length(rep_names)
    y_mpc2 = [y_mpc2; Fast_MPC.(char(rep_names(i))).output];  
end

load('e_2.mat');
y_mpc3 = [];
for i=1:1:length(rep_names)
    y_mpc3 = [y_mpc3; Fast_MPC.(char(rep_names(i))).output];  
end

load('e_3.mat');
y_mpc4 = [];
for i=1:1:length(rep_names)
    y_mpc4 = [y_mpc4; Fast_MPC.(char(rep_names(i))).output];   
end

load('e_4.mat');
y_mpc5 = [];
for i=1:1:length(rep_names)
    y_mpc5 = [y_mpc5; Fast_MPC.(char(rep_names(i))).output];  
end

load('e_5.mat');
y_mpc6 = [];
for i=1:1:length(rep_names)
    y_mpc6 = [y_mpc6; Fast_MPC.(char(rep_names(i))).output];  
end

load('0.mat');
y_mpc7 = [];
for i=1:1:length(rep_names)
    y_mpc7 = [y_mpc7; Fast_MPC.(char(rep_names(i))).output];  
end

figure(1)
boxes=[y_nom, y_mpc1 ,y_mpc2, y_mpc3, y_mpc4, y_mpc5, y_mpc6, y_mpc7] ;
%boxplot(boxes);
violin(boxes,'facecolor',[b;o;y;p;g;lightb;r;b2]);
hold on;
xL = get(gca,'XLim');
line(xL,[ysp ysp],'Color',g2,'Linewidth',2,'LineStyle','--');
ylabel('Final viable cell conc. cells/L')
xticklabels({'Nominal','M=1','M=1e-1','M=1e-2','M=1e-3','M=1e-4','M=1e-5','M=0'})
%legend('ODE generated','ODE + measure noise','Soft sensor')
hl = legend;
hl.String{3} = 'Set-point';
set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(gcf, 'Color', 'w');
savefig('cells_const')

load('VC.mat');
y_vc = [];
for i=1:1:length(rep_names)
    y_vc = [y_vc; Fast_MPC.(char(rep_names(i))).output];  
end

load('linear.mat');
y_lin = [];
for i=1:1:length(rep_names)
    y_lin = [y_lin; Fast_MPC.(char(rep_names(i))).output];   
end

load('exp.mat');
y_exp = [];
for i=1:1:length(rep_names)
    y_exp = [y_exp; Fast_MPC.(char(rep_names(i))).output];   
end

load('DE.mat');
y_de = [];
for i=1:1:length(rep_names)
    y_de = [y_de; Fast_MPC.(char(rep_names(i))).output];   
end


load('BSCI.mat');
y_bci = [];
for i=1:1:length(rep_names)
    y_bci = [y_bci; Fast_MPC.(char(rep_names(i))).output];   
end


load('BSCI2.mat');
y_bci2 = [];
for i=1:1:length(rep_names)
    y_bci2 = [y_bci2; Fast_MPC.(char(rep_names(i))).output];   
end


figure(2)
boxes=[y_nom ,y_vc, y_lin, y_exp, y_de, y_bci, y_bci] ;
%boxplot(boxes);
violin(boxes,'facecolor',[b;o;y;g;p;lightb;r]);
hold on;
xL = get(gca,'XLim');
line(xL,[ysp ysp],'Color',g2,'Linewidth',2,'LineStyle','--'); 
ylabel('Final viable cell conc. cells/L')
xticklabels({'Nominal','VC','Linear','Exponential','DE','BCI1','BCI2'})
%legend('ODE generated','ODE + measure noise','Soft sensor')
hl = legend;
hl.String{3} = 'Set-point';
set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(gcf, 'Color', 'w');
savefig('cells_comp')



 figure(3)
 plot(y_nom(1:80),'Color',b,'LineStyle',':','Marker','s','linewidth',1.5,'MarkerSize',8);
 hold on;
 %plot(y_mpc4(1:80),'blue--o','linewidth',1.5);
 plot(y_exp(1:80),'Color',p,'LineStyle',':','Marker','x','linewidth',1.5,'MarkerSize',8);
 plot(y_bci(1:80),'Color',r,'LineStyle',':','Marker','d','linewidth',1.5,'MarkerSize',8);
 plot(y_bci2(1:80),'Color',y,'LineStyle',':','Marker','o','linewidth',1.5,'MarkerSize',8);
 set(get(gca,'YLabel'),'String','Final viable cell conc. cells/L','fontsize',20)
 set(get(gca,'XLabel'),'String','Batch number','fontsize',20)
 xl=[0 80];
 yl=[2.84e9 2.84e9];
 line(xl,yl,'Color',g2,'LineStyle','--','linewidth',1.5);
 set(findall(gcf,'-property','FontSize'),'FontSize',24)
 legend('Nominal','Exponential','BCI1','BCI2','Set-point', 'FontSize',20,'Location','best');
 set(gcf, 'Color', 'w');
 

 
figure(5)
boxes=[y_nom ,y_mpc4 , y_bci, y_bci2] ;
%boxplot(boxes);
violin(boxes,'facecolor',[b;o;y;g]);
hold on;
xL = get(gca,'XLim');
line(xL,[ysp ysp],'Color',g2,'Linewidth',2,'LineStyle','--'); 
ylabel('Final viable cell conc. cells/L')
xticklabels({'Nominal','M=1e-3','BCI','BCI2'})
%legend('ODE generated','ODE + measure noise','Soft sensor')
hl = legend;
hl.String{3} = 'Set-point';
set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(gcf, 'Color', 'w');

sum(abs(y_mpc1-2.84e9))/length(y_mpc1)





