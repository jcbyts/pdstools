function [thresh, x, y] = interpThresh(inference, level)
% use linear interpolation to recover threshold from a psignifit inference
% [thresh, fitX, fitY] = interpThresh(inference, level)
% INPUTS:
% 	inference [struct] - output of BootstrapInference or BayesianInference from psignifit3.x
% 		level [double] - threshold level (default = .75)
% OUTPUTS:
% 	   thresh [double] - stimulus value that reaches [level] performance
% 	     fitX  [m x 1] - vector of stimulus values
% 		 fitY  [m x 1] - vector of PMF values (for plotting)



% 20140825 jly 	wrote it
% 20140907 jly  use psignifit pmf diagnostics so there are no 'core'
%               requirements. Uses only unique values so the function is
%               always monotonic.

import pdsa.* % matters if running from pdstools

thresh = [];

if nargin < 2 || isempty(level)
	level = .75;
	if nargin < 1
		help interpThresh
		return
	end
end

% switch inference.nafc
%     case 2
%         x = linspace(0,1,100);
%     case 1
% NOT HAPPY WITH THIS ... needs work
        x = linspace(min(inference.data(:,1)), max(inference.data(:,1))*10, 1000);
% end


% % Diagnostics of the point estimate
% if inference.gammaislambda
%     diag = Diagnostics ( inference.data, inference.params_estimate, ...
%         'sigmoid', inference.sigmoid, 'core', inference.core, 'nafc', inference.nafc, 'cuts', inference.cuts, 'gammaislambda' );
% else
%     diag = Diagnostics ( inference.data, inference.params_estimate, ...
%         'sigmoid', inference.sigmoid, 'core', inference.core, 'nafc', inference.nafc, 'cuts', inference.cuts );
% end
y = evaluate(x, inference);
% x = diag.pmf(:,1);
% y = diag.pmf(:,2);

% if things aren't monotonic because it saturates too early
if numel(unique(y)) < numel(y)
    uy = unique(y);
    nu = numel(unique(y));
    yInd = zeros(nu,1);
    for ii = 1:nu
        [~,yInd(ii)] = find(y==uy(ii), 1, 'first');
    end
    y = y(yInd);
    x = x(yInd);
end
if max(y) < level
    thresh = inf;
else
    thresh = interp1(y, x, level);
end


% old version with core ab requirements
% requirements:
% 	inference.core 	  MUST be set to 'ab' for this to work
%   inference.sigmoid MUST be 'logistic' or 'gauss'

% switch inference.sigmoid
% 	case 'logistic'
% 		sigmoid	= 'logistic';
% 	case 'gauss'
% 		sigmoid = 'norm';
% 	otherwise
% 		error('sigmoid must be logistic or gaussian')
% end
% 
% assert(strcmp(inference.core, 'ab'), 'run psignifit with "core" set to "ab"')
% 
% x = linspace(0, 1, 1e3);
% y = cdf(sigmoid, x, inference.params_estimate(1), inference.params_estimate(2));
% 
% switch inference.nafc
% 	case 1
% 		y = inference.params_estimate(4) + (1 - inference.params_estimate(3) - inference.params_estimate(4))*y;
% 	case 2
% 		y = .5 + (.5 - inference.params_estimate(3))*y;
% end
