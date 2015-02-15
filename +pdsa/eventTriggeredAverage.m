function [stAn, stAn_SD, widx] = eventTriggeredAverage(an, st, win)
% calculate a triggered average
% [stAn, stAn_SD] = eventTriggeredAverage(an, ev, win)
% INPUTS
%   an [n x m] analog data (number of samples by number of channels)
%   st [k x 1] vector of event indices to trigger on
%  win [1 x 2] left and right window edges (eg. [-500 500]) 
% OUTPUTS
%       stAn - triggered average
%    stAn_SD - triggered standard deviation

%%  Spike Triggered LFP
if ~exist('win', 'var')
    win = [-500 500]; % -500 ms to +500 ms
end

widx = win(1):win(2); % -500 ms to +500 ms
windowSize = numel(widx);

stAn    = zeros(windowSize, size(an, 2));
stAn_SD = zeros(windowSize, size(an, 2));

st(isnan(st)) = []; % remove impossible spikes
st(st > size(an,1)) = [];

widxs = bsxfun(@plus, st, widx);
widxs(widxs(:,1) <= 0,:) = []; % ignore very early spikes
widxs(widxs(:,end) > size(an,1),:) = []; % ignore very late spikes

for kCh = 1:size(an,2)
    l = an(widxs, kCh);
    l = reshape(l, [], windowSize);
    stAn(:, kCh) = mean(l);
    stAn_SD(:, kCh) = std(l)/ sqrt(size(l,1));
end
