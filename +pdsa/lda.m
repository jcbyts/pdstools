classdef lda < handle
    % LDA performs Fisher linear discriminant analysis
    % f = r'm + k 
    %
    % l = lda(X,Y)
    % l.train
    % Yhat=l.predict(Xnew)
    
    properties
        m@double
        k@double
    end
    
    properties(Access=private)
        mu1@double
        mu2@double
        C@double
        classes
        X
        Y
    end
    
    methods
        function l=lda(X,Y) % constructor
            
            
            l.classes=unique(Y);
            assert(numel(l.classes)==2, 'There are more than two classes. Fisher LDA will not work here')
            
            l.X = X;
            l.Y = Y;
            
            Y=Y==max(l.classes); % make it logical aligned to max class
            
            % get the required values
%             l.C=X'*X/numel(Y);
            l.C = cov(X,1);
            l.mu1=mean(X(Y,:))';
            l.mu2=mean(X(~Y,:))';

        end
        
        function train(l, rho)
            if nargin<2
                rho=1;
            end
            Chat=l.C*rho + (1-rho)*eye(size(l.C,1));
            Cinv = pinv(Chat);
            l.m = Cinv*(l.mu1-l.mu2);
            l.k = .5 * (l.mu2'*(Cinv*l.mu2)-l.mu1'*(Cinv*l.mu1));
%             l.m = Chat\(l.mu1-l.mu2);
%             l.k = .5 * (l.mu2'*(Chat\l.mu2)-l.mu1'*(Chat\l.mu1));
        end
        
        function trainCV(l, rhoSteps, Kfolds)
            rhos = linspace(0,1,rhoSteps);
            
            folds = cvpartition(numel(l.Y), 'Kfold', Kfolds);
            pc = nan(rhoSteps, Kfolds);
            I = eye(size(l.C,1));
            
            for r = 1:rhoSteps
                rho = rhos(r);
                
                for i = 1:Kfolds
                    
                    xtrain = l.X(folds.training(i),:);
                    ytrain  = l.Y(folds.training(i)) == max(l.Y);
                    
                    xtest = l.X(folds.test(i),:);
                    ytest  = l.Y(folds.test(i)) == max(l.Y);
                    
                    C_ = cov(xtrain,1);
                    Cdiag = diag(diag(C_));
                    mu1_ = mean(xtrain(ytrain,:))';
                    mu2_ = mean(xtrain(~ytrain,:))';
                    
                    Chat = C_*rho + (1-rho)*Cdiag;
                    Cinv = pinv(Chat);
                    
                    mhat = Chat\(mu1_-mu2_);
                    khat = .5 * (mu2_'*(Cinv*mu2_)-mu1_'*(Cinv*mu1_));
                    
                    pc(r,i) = mean(ytest == ((xtest*mhat + khat) > 0));
                end
            end
            
            [~, ind] = max(pc);
                   
            
        end
        
        function diagtrain(l, rho)
            if nargin<2
                rho=1;
            end
            Chat=l.C*rho + (1-rho)*diag(diag(l.C));
            l.m = Chat\(l.mu1-l.mu2);
            l.k = .5 * (l.mu2'*(Chat\l.mu2)-l.mu1'*(Chat\l.mu1));
        end
        
        function [mOut, kOut]=coefficients(l)
            mOut=l.m;
            kOut=l.k;
        end
        
        function [yhat, yclass] = predict(l, X)
           yhat = X*l.m + l.k;
           tmp = double(yhat > 0) + 1;
           yclass = l.classes(tmp);
        end
    end
end