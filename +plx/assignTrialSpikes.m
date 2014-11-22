function spikes = assignTrialSpikes(spikes, start, stop)
% Assign trial to spikes (depreciated)
% spikes = assignTrialSpikes(spikes, start, stop)
% inputs:
%   spikes [struct] - output from pls_getSpikes
%    start  [n x 1] - vector of trial start times
%     stop  [n x 1] - vector of trial stop times
% output:
%   spikes [struct] - modified to include trial subfield

% 20140601  jly     wrote it

ntrials = numel(start);
assert(ntrials==numel(stop), 'start and stop need the same number of trials')

spikes.trial = zeros(numel(spikes.time),1);
for ii = 1:ntrials
    idx = spikes.time > start(ii) & spikes.time < stop(ii);
    spikes.trial(idx) = ii;
end
