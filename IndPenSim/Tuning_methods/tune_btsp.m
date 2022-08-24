function [M,uc] = tune_btsp(mina, maxa, u_length, num_u, num_x, X2, Y2, ustart, invar, BETA, lv)
%Function to select the parameters of the diagonal suppresion matrix
%using confidence intervals in bootstrap calculations
%Copyright Carlos Alberto Duran Villalobos, February 2022-University of Manchester
lenx=size(X2,2);
n = u_length;

nb=10000;
[x,bootsamp] = bootstrp(nb,@(bootr)pls(bootr(:,1:lenx),bootr(:,lenx+1),lv),[X2,Y2]);
se=0;
for i=1:1:nb
    se=se+(BETA'-x(i,:)).^2;
end
see1=sqrt(se/(nb-1));
see=see1*2;  %2 times std.
Bbt=abs(BETA)'-see;

mvt_ind=[];
uc=[];
for i=1:num_u
    for k=ustart:n
        index=(k-1)*num_x + (k-ustart)*num_u + i + num_x+invar;
        mvt_ind=[mvt_ind index];
    end
end

M=ones(n*num_u,1);
maxBbt = max(Bbt(mvt_ind));
minBbt = min(Bbt(mvt_ind));
M_c = (Bbt(mvt_ind)-maxBbt)/(minBbt-maxBbt)*(maxa-mina);%controlled MVTs

for i=1:num_u
    M(ustart + (i-1)*n:n*i)= M_c((i-1)*(n-ustart+1)+1:i*(n-ustart+1));
end

% figure(1)
% yyaxis left
% bar(Bbt(mvt_ind))
% hold on;
% ylabel('MVT regression coefficient')
% xlabel('Time(h)')
% set(findall(gcf,'-property','FontSize'),'FontSize',24)
% yyaxis right
% plot([M(ustart:length(M)/num_u);M(ustart+length(M)/num_u:end)],'linewidth',2)
% ylabel('Suppresion matrix weight')
% set(findall(gcf,'-property','FontSize'),'FontSize',24)
end

