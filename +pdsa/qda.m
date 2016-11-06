classdef qda < handle
    % QDA implements Quadratic Discriminant Analysis
    % Solves the classification problem
    % f=r'Qr + r'm + k
    %
    % Call q=qda(X,Y);
    % q.train
    % q.predict(X)
    % q.coefficients
    properties
        Q@double
        m@double
        k@double
    end
    
    properties (Access=private)
        C1@double
        C2@double
        mu1@double
        mu2@double
        C1inv@double
        C2inv@double
    end
    
    methods
        
        function q=qda(X,Y)
            % QDA constructor
            
            classes=unique(Y);
            assert(numel(classes)==2, 'There are more than two classes. Fisher QDA will not work here')
            
            Y=Y==max(classes); % make it logical aligned to the biggest class ID
            
            q.C1=X(Y,:)'*X(Y,:)/sum(Y);
            q.C2=X(~Y,:)'*X(~Y,:)/sum(~Y);

            q.mu1=mean(X(Y,:))';
            q.mu2=mean(X(~Y,:))';
            
            q.C1inv=pinv(q.C1);
            q.C2inv=pinv(q.C2);
            
        end
        
        
        function train(q)
            q.Q = -.5*(q.C1inv-q.C2inv);
            c1m1=q.C1inv*q.mu1;
            c2m2=q.C2inv*q.mu2;
            q.k = -.5*(q.logdet(q.C1) - q.logdet(q.C2) - q.mu1'*c1m1 + q.mu2'*c2m2);
            q.m = c2m2 - c1m1;
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