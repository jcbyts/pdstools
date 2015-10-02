function [cnt,dt] = countSpikes(sptimes, ev1, ev2)
% Count spikes between two events
% v = countSpikes(sptimes, ev1, ev2)
% input
%       sptimes    [n x 1] - vector of spike times
%       ev1        [m x 1] - vector of event times to start count
%       ev2 [(p or 1) x 1] - vector of event times to stop count or
%                            scalar window size to count
% output
%       cnt - vector of spike counts
%       dt  - vector of time windows

% 20140416 jly wrote it

nE1 = numel(ev1);
nE2 = numel(ev2);

assert(nE2==1 | nE2==nE1, 'event 2 must have same # as event 1, or be a window size')

if nE2 == 1 && nE1~=nE2
    ev2 = ev1+ev2;
end

invalidEvents=isnan(ev1)|isnan(ev2);

cnt=nan(nE1,1);
dt=nan(nE1,1);

for ii = 1:nE1
    if invalidEvents(ii)
        continue
    end
    
    cnt(ii)=sum(sptimes>ev1(ii) & sptimes<ev2(ii));
    dt(ii)=ev2(ii)-ev1(ii);
end