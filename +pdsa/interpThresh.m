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
% requirements:
% 	inference.core 	  MUST be set to 'ab' for this to work
%   inference.sigmoid MUST be 'logistic' or 'gauss'

% 20140825 jly 	wrote it

import pdsa.* % matters if running from pdstools

thresh = [];

if nargin < 2 || isempty(level)
	level = .75;
	if nargin < 1
		help interpThresh
		return
	end
end

switch inference.sigmoid
	case 'logistic'
		sigmoid	= 'logistic';
	case 'gauss'
		sigmoid = 'norm';
	otherwise
		error('sigmoid must be logistic or gaussian')
end

assert(strcmp(inference.core, 'ab'), 'run psignifit with "core" set to "ab"')

x = linspace(0, 1, 1e3);
y = cdf(sigmoid, x, inference.params_estimate(1), inference.params_estimate(2));

switch inference.nafc
	case 1
		y = inference.params_estimate(4) + (1 - inference.params_estimate(3) - inference.params_estimate(4))*y;
	case 2
		y = .5 + (.5 - inference.params_estimate(3))*y;
end

thresh = interp1(y, x, level);
