
%% simulate some gaussian responses with different means and same covariance
n=20;

mu1=randn(1,n)+10;
mu2=randn(1,n)+8;

C=toeplitz([1 .5 zeros(1,n-2)])*10;

nTrials=200;
Y=rand(nTrials,1)<.5; % conditions
X=zeros(nTrials,n);
X(Y,:)=mvnrnd(mu1, C, sum(Y));
X(~Y,:)=mvnrnd(mu2, C, sum(~Y));


l=pdsa.lda(X,Y);
l.train
figure(1); clf
subplot(121)
histogram(X(Y,1)); hold on
histogram(X(~Y,1));
subplot(122)
yhat=l.predict(X);
histogram(yhat(Y)); hold on
histogram(yhat(~Y));

%% simulate responses to two classes with different covariances (QDA)
n=20;

mu1=zeros(1,n);
mu2=zeros(1,n);

C1=toeplitz([1 .5 .2 zeros(1,n-3)])*2;
C2=toeplitz([2 1 .2 zeros(1,n-3)])*2;
C1=C1'*C1;
C2=C2'*C2;

nTrials=200;
Y=rand(nTrials,1)<.5; % conditions
X=zeros(nTrials,n);
X(Y,:)=mvnrnd(mu1, C1, sum(Y));
X(~Y,:)=mvnrnd(mu2, C2, sum(~Y));


l=pdsa.lda(X,Y);
l.train
figure(1); clf
subplot(131)
histogram(X(Y,1)); hold on
histogram(X(~Y,1));
title('Single Response', 'FontWeight', 'Normal')
subplot(132)
yhat=l.predict(X);
histogram(yhat(Y)); hold on
histogram(yhat(~Y));
title('LDA', 'FontWeight', 'Normal')
q=pdsa.qda(X,Y);
q.train;
[Q,m,k]=q.coefficients;
subplot(133)
Yhat=q.predict(X);
histogram(Yhat(Y)); hold on
histogram(Yhat(~Y));
title('QDA', 'FontWeight', 'Normal')