
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


[m,k]=pdsa.lda(X,Y);

figure(1); clf
subplot(121)
histogram(X(Y,1)); hold on
histogram(X(~Y,1));
subplot(122)
histogram(X(Y,:)*m+k); hold on
histogram(X(~Y,:)*m+k);

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


[m,k]=pdsa.lda(X,Y);

figure(1); clf
subplot(131)
histogram(X(Y,1)); hold on
histogram(X(~Y,1));
subplot(132)
histogram(X(Y,:)*m+k); hold on
histogram(X(~Y,:)*m+k);

q=pdsa.qda(X,Y);
q.train;
[Q,m,k]=q.coeficients;
subplot(133)
Yhat=q.predict(X);
histogram(Yhat(Y)); hold on
histogram(Yhat(~Y));