function MCCV=plsmccv(X,y,A,N,ratio,OPT)
%+++ Monte Carlo Cross Validation for PLS regression.
%+++ Input:  X: m x n  (Sample matrix)
%            y: m x 1  (measured property)
%            A: The maximal number of latent variables for
%            cross-validation
%            N: The number of Monte Carlo Simulation.
%        ratio: The ratio of calibration samples to the total samples.
%               pareto,minmax,center or none.
%          OPT: =1 : print process.
%               =0 : don't print process.
%+++ Output: Structural data: CV
%+++ Ref: Q.S. Xu, Y.Z. Liang, 2001.Chemo Lab,1-11
%+++ Edited by: Carlos Duran- Manchester University 2017

function [B,Wstar,T,P,Q,W,R2X,R2Y]=plsnipals(X,Y,A)
%+++ The NIPALS algorithm for both PLS-1 (a single y) and PLS-2 (multiple Y)
%+++ X: n x p matrix
%+++ Y: n x m matrix
%+++ A: number of latent variables



varX=sum(sum(X.^2));
varY=sum(sum(Y.^2));
for i=1:A
    error=1;
    u=Y(:,1);
    niter=0;
    while (error>1e-8 && niter<1000)  % for convergence test
        w=X'*u/(u'*u);
        w=w/norm(w);
        t=X*w;
        q=Y'*t/(t'*t);  % regress Y against t;
        u1=Y*q/(q'*q);
        error=norm(u1-u)/norm(u);
        u=u1;
        niter=niter+1;
    end
    p=X'*t/(t'*t);
    X=X-t*p';
    Y=Y-t*q';
    
    %+++ store
    W(:,i)=w;
    T(:,i)=t;
    P(:,i)=p;
    Q(:,i)=q;
    
end

%+++ calculate explained variance
R2X=diag(T'*T*P'*P)/varX;
R2Y=diag(T'*T*Q'*Q)/varY;

Wstar=W*(P'*W)^(-1); 
B=Wstar*Q';
Q=Q';

end

if nargin<7;OPT=0;end
if nargin<6;ratio=0.8;end
if nargin<5;N=1000;end
if nargin<3;A=2;end


tic;
[Mx,Nx]=size(X);
A=min([size(X) A]);
nc=floor(Mx*ratio);
nv=Mx-nc;
yytest=[];yycal=[];
YR=[];YC=[];Trace=[];
TrainIndex=[];
TestIndex=[];
for i=1:N
    
    % generate sub-training set and sub-test set
    index=randperm(Mx);
    calk=index(1:nc);
    testk=index(nc+1:Mx);
    Xcal=X(calk,:);ycal=y(calk);
    Xtest=X(testk,:);ytest=y(testk);
    
    %   data pretreatment
    [Xs,xpara1,xpara2]=zscore(Xcal);
    xpara2(xpara2<=1e-12)=1e-12;
    [ys,ypara1,ypara2]=zscore(ycal);   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[P,Q,T,U,BETA,PCTVAR,MSE,stats]=plsregress(Xs,ys,A);
    %B=BETA(2:end,:);
    %Wstar=stats.W;
    [B,Wstar,T,P,Q]=plsnipals(Xs,ys,A);
 
    yp=[];yc=[];
    for j=1:A
        B=Wstar(:,1:j)*Q(1:j);
        %+++ calculate the coefficient linking Xcal and ycal.
        C=ypara2*B./xpara2';
        coef=[C;ypara1-xpara1*C;];
        %+++ predict
        Xteste=[Xtest ones(size(Xtest,1),1)];
        Xcale=[Xcal ones(size(Xcal,1),1)];
        ycal_p=Xcale*coef;       ytest_p=Xteste*coef;
        yp=[yp ytest_p];         yc=[yc ycal_p];     
    end
     
    YR=[YR;yp];yytest=[yytest;ytest];
    
    e1m=sum((yc-repmat(ycal,1,A)).^2);
    e2m=sum((yp-repmat(ytest,1,A)).^2);
    
    e1=sqrt(e1m/nc);
    e2=sqrt(e2m/nv);
    
    Trace=[Trace;[e1 e2]];
    
    if OPT==1;if mod(i,100)==0;fprintf('The %dth sampling for MCCV finished.\n',i);end;end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_test=YR-repmat(yytest,1,A);
ntotal=length(yytest);
error2=error_test.^2;
error2_MEAN=sum(error2)/ntotal;
error2_SD= sqrt(sum((error2-repmat(mean(error2),ntotal,1)).^2)/(ntotal-1)); % unbiased estimator


PRESS=sum(error2);
cv=sqrt(PRESS/N/nv);
[RMSECV_min,index]=min(cv);
index=min(index);
SST=sumsqr(yytest-mean(yytest));
for i=1:A
  SSE=sumsqr(YR(:,i)-yytest);
  Q2(i)=1-SSE/SST;
end

indexSD=find(error2_MEAN <= min(error2_MEAN)+error2_SD(index));
indexSD=min(indexSD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%+++ output  %%%%%%%%%%%%%%%%
MCCV.MC_para=[N ratio];
MCCV.Ypred=YR;
MCCV.Ytrue=yytest;
MCCV.RMSECV=cv;
MCCV.RMSECV_min=RMSECV_min;
MCCV.Q2=Q2;
MCCV.Q2_max=Q2(index);
MCCV.optLV=index;
MCCV.RMSEF=Trace(:,1:A);
MCCV.RMSEP=Trace(:,A+1:end);
MCCV.note='****** The following is based on global min MSE + 1SD';
MCCV.RMSECV_min_1SD=cv(indexSD);
MCCV.Q2_max_1SD=Q2(indexSD);
MCCV.optLV_1SD=indexSD;
MCCV.time=toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


