% To organize PLS matrices based on the location of the control point
% Author: Carlos Duran - Manchester University 2017
% Last modified: March , 2021

function [Xtrn,xp,Pn,Wn,En,xmeann,xstdn]=orgplsmpc(X2,P,W,E,xmean,xstd,num_x,num_u,bt,ustart,xc,cp,invar)


Xtrn=X2(:,1:invar);
xmeann=xmean(1,1:invar);
xstdn=xstd(1,1:invar);
Pn=P(1:invar,:);
Wn=W(1:invar,:);
En=E(:,1:invar);
xp=xc(1:invar);

%xpast before ustart
for i=1:num_x
    for k=1:ustart-1
        index=(k-1)*num_x+i+invar;
        xp=[xp xc(index)];
        Xtrn=[Xtrn X2(:,index)];
        xmeann=[xmeann xmean(1,index)];
        xstdn=[xstdn xstd(1,index)];
        Pn=[Pn;P(index,:)];
        Wn=[Wn;W(index,:)];
        En=[En,E(:,index)];
    end
end
%xpast before cp
for i=1:num_x
    for k=ustart:cp-1
        index=(k-1)*num_x +(k-ustart)*num_u + i+invar;
        xp=[xp xc(index)];
        Xtrn=[Xtrn X2(:,index)];
        xmeann=[xmeann xmean(1,index)];
        xstdn=[xstdn xstd(1,index)];
        Pn=[Pn;P(index,:)];
        Wn=[Wn;W(index,:)];
        En=[En,E(:,index)];
    end
end
%upast before cp
for i=1:num_u
    for k=ustart:cp-1
        index=(k-1)*num_x + (k-ustart)*num_u + i + num_x+invar;
        xp=[xp xc(index)];
        Xtrn=[Xtrn X2(:,index)];
        xmeann=[xmeann xmean(1,index)];
        xstdn=[xstdn xstd(1,index)];
        Pn=[Pn;P(index,:)];
        Wn=[Wn;W(index,:)];
        En=[En,E(:,index)];
    end
end
%ufuture after cp
for i=1:num_u
    for k=cp:bt
        index=(k-1)*num_x + (k-ustart)*num_u + i + num_x+invar;
        Xtrn=[Xtrn X2(:,index)];
        xmeann=[xmeann xmean(1,index)];
        xstdn=[xstdn xstd(1,index)];
        Pn=[Pn;P(index,:)];
        Wn=[Wn;W(index,:)];
        En=[En,E(:,index)];
    end
end

%xfuture
for i=1:num_x
    for k=cp:bt
        index=(k-1)*num_x +(k-ustart)*num_u + i+invar;
        Xtrn=[Xtrn X2(:,index)];
        xmeann=[xmeann xmean(1,index)];
        xstdn=[xstdn xstd(1,index)];
        Pn=[Pn;P(index,:)];
        Wn=[Wn;W(index,:)];
        En=[En,E(:,index)];
    end
end
end


