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
    end
    
    methods
        function l=lda(X,Y,varargin) % constructor
            
            
            classes=unique(Y);
            assert(numel(classes)==2, 'There are more than two classes. Fisher LDA will not work here')
            
            Y=Y==max(classes); % make it logical aligned to max class
            
            % get the required values
            l.C=X'*X/numel(Y);
            l.mu1=mean(X(Y,:))';
            l.mu2=mean(X(~Y,:))';

        end
        
        function train(l)
            l.m = l.C\(l.mu1-l.mu2);
            l.k = .5 * (l.mu2'*(l.C\l.mu2)-l.mu1'*(l.C\l.mu1));
        end
        
        function [mOut, kOut]=coefficients(l)
            mOut=l.m;
            kOut=l.k;
        end
        
        function yhat=predict(l, X)
           yhat=X*l.m + l.k;
        end
    end
end