classdef qda < handle
    % QDA implements Quadratic Discriminant Analysis
    % Solves the classification problem
    % f=r'Qr + r'm + k
    %
    % Call q=qda(X,Y);
    % q.train(1,1) - QDA
    % q.train(1,0) - LDA
    % q.train(0,1) - QDA Identity covariance
    % q.train(.5,1) - QDA ridge regularization
    % q.train(.5,.5) - ridge regularization and mixture between
    % q.predict(X)
    % q.coefficients
    properties
        Q@double
        m@double
        k@double
    end
    
    properties (Access=private)
        C@double
        C1@double
        C2@double
        mu1@double
        mu2@double
    end
    
    methods
        
        function q=qda(X,Y)
            % QDA constructor
            
            classes=unique(Y);
            assert(numel(classes)==2, 'There are more than two classes. Fisher QDA will not work here')
            
            Y=Y==max(classes); % make it logical aligned to the biggest class ID
            
            q.C1=X(Y,:)'*X(Y,:)/sum(Y);
            q.C2=X(~Y,:)'*X(~Y,:)/sum(~Y);
            q.C=X'*X/numel(Y);
            
            q.mu1=mean(X(Y,:))';
            q.mu2=mean(X(~Y,:))';
                        
        end
        
        
        function train(q, gamma, alpha, r, diagonalize)
            % Train QDA
            % train(q, gamma, alpha, r, diagonalize)
            %
            % regularize in two ways
            % 1) take steps between QDA and LDA
            %   Chat1 = alpha*Chat1 + (1-alpha)*Chat
            %
            %   if alpha=1, QDA
            %   if alpha=0, LDA
            %
            % 2) use ridge on the covariance
            %   Chat = gamma*Chat + (1-gamma)*I
            %
            %   if gamma=0, Chat = Identity
            
            if nargin<5 || isempty(diagonalize)
                diagonalize=0;
            end
            
            if nargin<4 || isempty(r)
                r=max(size(q.C));
            end
            
            if nargin<3 || isempty(alpha)
                alpha=1;
            end
            
            if nargin<2 || isempty(gamma)
                gamma=1;
            end
            
            n=size(q.C,1);
            
            if diagonalize
                Chat=diag(diag(q.C));
                C1hat=diag(diag(q.C1));
                C2hat=diag(diag(q.C2));
            else
                Chat=q.C;
                C1hat=q.C1;
                C2hat=q.C2;
            end
            
            % ridge regularization for estimation of the sample covariance
            Chat=gamma*Chat + (1-gamma)*speye(n);
            C1hat=gamma*C1hat + (1-gamma)*speye(n);
            C2hat=gamma*C2hat + (1-gamma)*speye(n);
            
            % step between QDA and LDA
            C1hat=alpha*C1hat + (1-alpha)*Chat;
            C2hat=alpha*C2hat + (1-alpha)*Chat;
            
            % low rank approximation of the inverse covariance
            if r<max(size(q.C))
                [u, s, v]=svd(C1hat);
                sd=diag(s);
                shat=diag(1./sd(1:r));
                C1invhat=u(:,1:r)*shat*v(:,1:r)';
                [u, s, v]=svd(C2hat);
                sd=diag(s);
                shat=diag(1./sd(1:r));
                C2invhat=u(:,1:r)*shat*v(:,1:r)';
            else
                C1invhat=pinv(C1hat);
                C2invhat=pinv(C2hat);
            end
            % Regularized estimates of Q, m, k
            q.Q = -.5*(C1invhat-C2invhat);
            
            c1m1=(C1invhat)*q.mu1;
            c2m2=(C2invhat)*q.mu2;
            
            q.m = c2m2 - c1m1;
                        
            ldc1=q.logdet(C1hat);
            ldc2=q.logdet(C2hat);
            
            q.k = -.5*(ldc1 - ldc2 - q.mu1'*c1m1 + q.mu2'*c2m2);
            
            %             q.Q = -.5*(q.C1inv-q.C2inv);
            %             c1m1=q.C1inv*q.mu1;
            %             c2m2=q.C2inv*q.mu2;
            %             q.k = -.5*(q.logdet(q.C1) - q.logdet(q.C2) - q.mu1'*c1m1 + q.mu2'*c2m2);
            %             q.m = c2m2 - c1m1;
        end
        
        
        function [qOut, mOut, kOut]=coefficients(q)
            qOut=q.Q;
            mOut=q.m;
            kOut=q.k;
        end
        
        function yhat=predict(q,X)
            n=size(X,1);
            yhat=nan(n,1);
            for i=1:n
                yhat(i)=X(i,:)*q.Q*X(i,:)' + X(i,:)*q.m + q.k;
            end
            
            if norm(q.Q)==0 % if the quadratic term is zero, flip sign
                yhat=-yhat;
            end
        end
        
    end
    
    methods (Static)
        
        function x = logdet(A)
            % LOGDET - computes the log-determinant of a matrix A
            %
            % x = logdet(A);
            %
            % This is faster and more stable than using log(det(A))
            %
            % Input:
            %     A NxN - A must be sqaure, positive semi-definite
            
            x = 2*sum(log(diag(chol(A))));
            
            
        end
    end
end