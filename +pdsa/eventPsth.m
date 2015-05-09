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


[ev, trid] = sort(ev(:));


be = win(1):bs:win(2); 
bc = be(1:end-1)+bs/2;
nbins = numel(bc);

% normalize smoothing kernel
if ~isempty(skern)
    k = numel(skern);
    be = (win(1)-bs*k):bs:(win(2)+bs*k);
    bc = be(1:end-1)+bs/2;
    nbins = numel(bc);
    goodindex = bc>=win(1) & bc <=win(2);
    skern = skern/sum(skern);
    znorm = filtfilt(skern,1,ones(nbins+1,1));
end

validEvents = ~isnan(ev);
nEvents = sum(validEvents);

% for loop version (super slow)
chunkSize = 1e3;
if nEvents > chunkSize
    evstart = (1:chunkSize:nEvents)';
    evend   = [chunkSize:chunkSize:nEvents nEvents]';
    
    tspcnt  = zeros(nEvents, nbins);
    for kk = 1:numel(evstart)
        idx         = evstart(kk):evend(kk);
        evtmp       = ev(validEvents(idx))';
        spindx      = sptimes > (evtmp(1) - win(1)) & sptimes < (evtmp(end) + win(2));
%         spindx      = 1:numel(sptimes);

        tmpspcnt    = bsxfun(@minus, sptimes(spindx), evtmp);
        tmpspcnt    = histc(tmpspcnt, be)';
        if any(tmpspcnt(:))
            tspcnt(idx,:)  = tmpspcnt(:,1:nbins);
%             tspcnt(trid(idx),:)  = tmpspcnt(:,1:nbins);
        end
    end
    
%     for ii = 1:nEvents
%         ievent = validEvents(ii);
%         st = sptimes - ev(ievent);
%         st = st(st>be(1) & st < be(end));
%         spcnt = histc(st, be);
%         smcnt = spcnt;
%         
%         if any(spcnt)
%             tspcnt(ii,:)  = smcnt(1:nbins);
%         end
%     end
else
    
    % bsxfun version
    tmpspcnt = bsxfun(@minus, sptimes, ev(validEvents)');
    
    tspcnt    = histc(tmpspcnt, be)';
    tspcnt(:,end) = [];
end


% if smoothing needs to be done
if ~isempty(skern)
%     tspcnt = filter(skern, 1, tspcnt, [], 2);
%     nz = filter(skern, 1, ones(size(tspcnt)), [], 2);

    tspcnt = filtfilt(skern, 1, tspcnt')';
    nz = filtfilt(skern, 1, ones(size(tspcnt))')';
    tspcnt = tspcnt./nz;
    tspcnt = tspcnt(:, goodindex);
    bc = bc(goodindex);
end

% get summary statistics
m = mean(tspcnt)/bs;
v = var(tspcnt)/bs;
s = std(tspcnt)/sqrt(nEvents)/bs;

tspcntOut = tspcnt;
% if nargout > 4
%     tspcntOut = zeros(size(tspcnt));
%     tspcntOut(trid,:) = tspcnt;
% end

    
