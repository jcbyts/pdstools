function [xc, lags, r] = xcorrShuffle(sp1, sp2, ev, ml, win, bs, shuf)
% Compute shuffle corrected cross-correlogram 
% [xc, lags, r] = xcorrShuffle(sp1, sp2, ev, ml, win, bs, shuf)
% inputs
%       sp1     [n x 1] - vector of spike times for neuron 1
%       sp2     [m x 1] - vector of spike times for neuron 2
%       ev      [k x 1] - vector of event times to align to
%       ml      [1 x 1] - scalar (max number of lags for ccg)
%       win     [1 x 2] - range to count spikes from after aligning
%       bs      [1 x 1] - bin size
%       shuf    [1 x 1] - number shuffle permutations (0 is no shuffle)

% 20140416 jly wrote it

if nargin < 7
    shuf = 1;
    if nargin < 6
        if nargin < 5
            win = [-mean(diff(sp1)) 10*mean(diff(sp1))];
            if nargin < 4
                ml = 200;
                if nargin < 3
                    help xcorrShuffle
                    return
                end
            end
            bs = diff(win)/100;
        end
    end
end

validEvents = find(~isnan(ev));

nEvents = numel(validEvents);
shuffleEvents = validEvents(randi(nEvents, [nEvents 1]));

be  = (-(ml*bs):bs:(ml*bs))-(bs/2);
lags  = be(1:end-1)+bs/2;
nbins = numel(lags);
xct    = zeros(nEvents, nbins);
tspcnt1 = zeros(nEvents,1);
tspcnt2 = zeros(nEvents,1);
tspcnt3 = zeros(nEvents,1);
tspcnt4 = zeros(nEvents,1);
xcshuffle = zeros(nEvents, nbins);
for ii = 1:nEvents
    ievent = validEvents(ii);
    st1 = sp1 - ev(ievent);
    st2 = sp2 - ev(ievent);
    
    st1 = st1(st1 >= win(1) & st1 <= win(2));
    st2 = st2(st2 >= win(1) & st2 <= win(2));
    
    tspcnt1(ii) = numel(st1);
    tspcnt2(ii) = numel(st2);
    
    % differences between each spike time
    d = bsxfun(@minus, repmat(st1(:), 1, numel(st2)), st2(:)');
    h = histc(d(:), be);
    xct(ii,:) = h(1:end-1);
    
    if shuf
        h = 0;
        tmpcnt3 = 0;
        tmpcnt4 = 0;
        for jj = 1:shuf
            % shuffled trials for neuron 1 and 2
            st3 = sp1 - ev(randsample(shuffleEvents,1)); % shuffled spikes
            st3 = st3(st3 >= win(1) & st3 <= win(2));
            st4 = sp2 - ev(randsample(shuffleEvents,1)); % shuffled spikes
            st4 = st4(st4 >= win(1) & st4 <= win(2));
            
            tmpcnt3 = tmpcnt3 + numel(st3);
            tmpcnt4 = tmpcnt4 + numel(st4);
            % differences between each spike time
            d = bsxfun(@minus, repmat(st1(:), 1, numel(st4)), st4(:)');
            a = histc(d(:), be); h = h(:) + a(:);
        end
        xcshuffle(ii,:) = h(1:end-1)/shuf;
        tspcnt3(ii) =  tmpcnt3/shuf;
        tspcnt4(ii) = tmpcnt4/shuf;
    end


end

xc = mean(xct);
pr = corrcoef(tspcnt1, tspcnt2);
r = pr(2);
if shuf
    xc = xc - mean(xcshuffle);
    ps = corrcoef(tspcnt3, tspcnt4);
%     [ps,p] = corrcoef(tspcnt1-tspcnt3, tspcnt2-tspcnt4);
    r = r - ps(2);
end

