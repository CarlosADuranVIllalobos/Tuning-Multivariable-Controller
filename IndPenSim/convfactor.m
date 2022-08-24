%% convfactor.m
% Function to convert struct indpensim_run.m for PLS

%% Copyright
% Carlos Alberto Duran Villalobos July 2020
% University of Manchester.

%% Code
function [x,y] = convfactor(xref,ustart)
x=[];
y= xref.P.y(end); %penicilin concentration for batches
steps=size(xref.DO2.y',2);
%batch wise
for k=1:5:steps
    x=[x,xref.DO2.y(k),xref.T.y(k),xref.pH.y(k),xref.V.y(k),xref.CO2outgas.y(k),xref.O2.y(k)];

    if(k>=(ustart-1)*5)
        x=[x,xref.Fs.y(k),xref.Fpaa.y(k)];
    end
    
end
end
 