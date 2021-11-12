%% convfactor.m
% Function to convert struct indpensim_run.m for PLS

%% Copyright
% Carlos Alberto Duran Villalobos July 2017
% University of Manchester.

%% Code
function [X,Y] = convfactor(Xref,ustart)
X=[];
% Y0(1); % total cell density
% Y0(2); % viable cells
% Y0(3); % dead cells
% Y0(4);  % glucose concentration
% Y0(5);  % glutamie concentration 
% Y0(6);  % lactate concentration
% Y0(7);  % ammonia concentration
% Y0(8); % unknown inhibitor concentration
% Y0(9);  % volume
% Y0(10); % glucose Flowrate set point

Y= Xref(end,2); %Final viable cell count
len=length(Xref);
 %batch wise
 for k=1:1:len
        X=[X, Xref(k,4), Xref(k,5), Xref(k,6), Xref(k,7), Xref(k,9)]; 
       % X=[X,Xref.DO2.y(k),Xref.T.y(k),Xref.pH.y(k),Xref.V.y(k),Xref.O2.y(k)]; 
        if(k>=(ustart))
        X=[X,Xref(k,10)];
        end     
 end
end
 