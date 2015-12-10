function eyeposF = upsampleEyepos(eyepos, Fs)
% up-sample and spline interpolate eye position to sampling rate
% eyeposF = upsampleEyepos(eyepos, Fs)
% eyepos is cell-array of [time x y] eye position vectors

if nargin < 2
    Fs = 1e3;
end

bs = 1/Fs;

nTrials = numel(eyepos);

eyeposF = cell(nTrials, 1);

for kTrial = 1:nTrials
    
    start = min(eyepos{kTrial}(:,1));
    stop  = max(eyepos{kTrial}(:,1));
    
    time_current = eyepos{kTrial}(:,1);
    
    time_new  = (start:bs:stop)';
    nsamples  = numel(time_new);
    x_new     = nan(nsamples, 1);
    y_new     = nan(nsamples, 1);
    
    
    dt = abs(bsxfun(@minus, time_new',time_current));
    
    [~, idx] = min(dt, [], 2);
    x_new(idx) = eyepos{kTrial}(:,2);
    y_new(idx) = eyepos{kTrial}(:,3);
    
    dx = mean([max(diff(find(diff(isnan(x_new))==1))) mean(diff(find(diff(isnan(x_new))==1)))]);
    dy = mean([max(diff(find(diff(isnan(x_new))==1))) mean(diff(find(diff(isnan(x_new))==1)))]);
    
    
    if isnan(dx) || isnan(dy)
        continue
    end
    
    dx = max(dx, 14);
    dy = max(dy, 14);
    % upsample and smooth
    x_filt = filtfilt(ones(ceil(dx),1)/ceil(dx), 1, repnan(x_new, 'spline'));
    y_filt = filtfilt(ones(ceil(dy),1)/ceil(dy), 1, repnan(y_new, 'spline'));
    
    eyeposF{kTrial} = [time_new(:) x_filt(:) y_filt(:)];
    
end

function [x] = repnan(x,method)
% REPNAN replaces NaN values in a 1D array.
% 
%% Syntax
% 
% x = repnan(x); 
% x = repnan(x,method); 
% 
%% Description 
% 
% x = repnan(x) returns x sans NaNs. 
% 
% x = repnan(x,method) specifies a method for replacing the original x's
% NaNs. Methods can be 'nearest', 'linear', 'spline', 'pchip', 'cubic',
% 'v5cubic', 'next', or 'previous'.  The 'next' options replaces NaNs with
% the next non-NaN value in x.  The 'previous' option replaces NaNs with
% the previous non-NaN value in x.  Default is 'linear'.  
% 
%% Author Info
% 
% Written by Chad A. Greene of the University of Texas at Austin's 
% Institute for Geophysics, October 31, 2014. Feel free to visit Chad 
% over at http://www.chadagreene.com. 
% 
% See also: interp1, find, NaN, and isnan.


%% Input check: 

assert(isvector(x)==1,'The repnan function requires that x must be a 1D vector.')

%% Use linear 1D interpolation by default: 

if nargin == 1
    method = 'linear'; 
end
    
%% Transpose to row vector if necessary: 

StartedColumn = false; 
if iscolumn(x) 
    x = x'; 
    StartedColumn = true; 
end

%% Perform interpolation: 

switch lower(method)
    case 'next'
        for k = find(isnan(x))
            try
                x(k) = x(find(1:length(x)>k & isfinite(x),1,'first'));
            end
        end
        
    case {'last','prev','previous'} 
        for k = find(isnan(x))
            try
                x(k) = x(find(1:length(x)<k & isfinite(x),1,'last'));
            end
        end
        
    case {'linear','cubic','nearest','spline','pchip','v5cubic'} 
            x(isnan(x)) = interp1(find(~isnan(x)),x(~isnan(x)),find(isnan(x)),method);
            
    otherwise
        error('Unrecognized method of interpolation.') 
end

%% Recolumnate x if user entered x as a column vector: 

if StartedColumn
    x = x'; 
end



