function [m,s,bc,v, tspcntOut] = eventPsth(sptimes, ev, win, bs, skern)
% Make Peri-Stimulus-Time-Histogram (PSTH) aligned to stimulus events
% [m, s, bc, v, tspcnt] = eventPsth(sptimes, ev, win, bs, skern)
%   inputs
%       sptimes [n x 1] - vector of spike times
%       ev      [m x 1] - vector of event times
%       win     [1 x 2] - window around events to count spikes in
%       bs      [1 x 1] - scalar bin size
%       skern   [1 x k] - smoothing kernel
%       works in any units as long as all variable have the same units.

% 20140415 jly wrote it

if nargin < 5
    skern = [];
    if nargin < 4
        if nargin < 3
            win = [-mean(diff(sptimes)) 10*mean(diff(sptimes))];
            if nargin < 2
                help eventPsth
                return
            end
        end
        bs = diff(win)/100;
    end
end

ev=ev(:);
binfun = @(t) (t == 0) + ceil(t/bs);

be = win(1):bs:win(2);
bc = be(1:end-1)+bs/2;

% normalize smoothing kernel
if ~isempty(skern)
    k = numel(skern);
    be = (win(1)-bs*k):bs:(win(2)+bs*k);
    bc = be(1:end-1)+bs/2;
    goodindex = bc>win(1) & bc <win(2);
    skern = skern/sum(skern);
end

validEvents = ~isnan(ev);
nTrials=numel(ev);

nEvents = sum(validEvents);

sbn = [];
str = [];
assert(nEvents<2e3, 'too many events for this to run fast!')
for kEvent=1:nTrials
    if ~validEvents(kEvent)
        continue
    end
    spo = sptimes(sptimes > ev(kEvent) + be(1) & sptimes < ev(kEvent) + be(end))- ev(kEvent);
    sbn = [sbn; binfun(spo- be(1))]; %#ok<AGROW>
    str = [str; ones(numel(spo),1)*kEvent]; %#ok<AGROW>
end
tspcnt = full(sparse(str, sbn, 1, nTrials, binfun(be(end)-be(1))));
tspcnt(~validEvents,:)=nan;

% if smoothing needs to be done
if ~isempty(skern)
    %     tspcnt = filter(skern, 1, tspcnt, [], 2);
    %     nz = filter(skern, 1, ones(size(tspcnt)), [], 2);    
    tspcnt(validEvents,:) = filtfilt(skern, 1, tspcnt(validEvents,:)')';
    nz = filtfilt(skern, 1, ones(size(tspcnt))')';
    tspcnt(validEvents,:) = tspcnt(validEvents,:)./nz(validEvents,:);
    tspcnt = tspcnt(:, goodindex);
    bc = bc(goodindex);
end

% get summary statistics
m = mean(tspcnt(validEvents,:))/bs;
v = var(tspcnt(validEvents,:))/bs;
s = std(tspcnt(validEvents,:))/sqrt(nEvents)/bs;

tspcntOut = tspcnt;