function [v, id, dt] = countSpikes(sptimes, ev1, ev2)
% Count spikes between two events
% v = countSpikes(sptimes, ev1, ev2)
% input
%       sptimes    [n x 1] - vector of spike times
%       ev1        [m x 1] - vector of event times to start count
%       ev2 [(p or 1) x 1] - vector of event times to stop count or
%                            scalar window size to count
% output
%       v  - vector of spike counts
%       id - index of ev1 that count comes from

% 20140416 jly wrote it


validEventStarts = find(~isnan(ev1));
validEventStops  = find(~isnan(ev2));
nE1 = numel(validEventStarts);
nE2 = numel(validEventStops);

if nE2 == 1 && nE1~=nE2
    e1 = ev1(validEventStarts);
    e2 = e1+ev2;
else
    e1   = ev1(validEventStarts);
    tmp  = ev2(validEventStops);
    e2   = nan(numel(e1),1);
    for ii = 1:nE1
        if ii == nE1
            if tmp(ii) > e1(ii)
                e2(ii) = tmp(ii);
            else
                e1(ii) = nan;
            end
            break
        end
        if ii == nE2
            break
        end
        eid = find(tmp>e1(ii) & tmp<e1(ii+1), 1);
        if ~isempty(eid)
            e2(ii) = tmp(eid);
        else
            e1(ii) = nan;
        end
    end
    
end

good = ~isnan(e1);
e1 = e1(good);
e2 = e2(good);
nE1 = sum(good);
dt = e2-e1;
id = validEventStarts(good);
v = zeros(nE1,1);
for ii = 1:nE1
   v(ii) = sum(sptimes > e1(ii) & sptimes < e2(ii));
end
