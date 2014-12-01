function [spikes, h] = stableSpikeRange(spikes, getInput)
% find range with stable recordings
% [spikes, h] = stableSpikeRange(spikes, getInput)

if nargin < 2
    getInput = false;
end
units = unique(spikes.id);
nUnits = numel(units);

h = nan(nUnits,1);
if isfield(spikes, 'goodRange')
    for ii = 1:nUnits
        stimes = spikes.time(spikes.id==units(ii));
        bs = 30; % seconds
        bins = 0:bs:max(stimes);
        spcnt = histc(stimes, bins)/bs;
        x = [bins(:) ones(numel(bins),1)];
        wls = (x'*x)\(x'*spcnt);
        rhat = x*wls;
        h(ii) = figure(ii); clf; hold all
        plot(bins, spcnt, 'o', bins, rhat, 'r')
        xlabel('time (seconds)')
        ylabel('spike rate')
        title(sprintf('unit: %d', units(ii)))
        idx = bins>spikes.goodRange(units(ii),1) & bins<spikes.goodRange(units(ii),2);
        plot(bins(idx), spcnt(idx), 'r.')
    end
   return 
end

spikes.goodRange = nan(numel(spikes.snr),2);

for ii = 1:nUnits
    stimes = spikes.time(spikes.id==units(ii));
    bs = 30; % seconds
    bins = 0:bs:max(stimes);
    spcnt = histc(stimes, bins)/bs;
    x = [bins(:) ones(numel(bins),1)];
    wls = (x'*x)\(x'*spcnt);
    rhat = x*wls;
    resids = spcnt - rhat; 
    h(ii) = figure(ii); clf; hold all
	plot(bins, spcnt, 'o', bins, rhat, 'r')

    aboveThreshold = resids<-(std(spcnt)/2);
    spanLocs = bwlabel(aboveThreshold);   %identify contiguous ones
    spanLength = regionprops(spanLocs, 'area');  %#ok<MRPBW> %length of each span
    spanLength = [ spanLength.Area];
    goodSpans = find(spanLength>=10);   %get on ly spans of 5+ points
    allInSpans = (ismember(spanLocs, goodSpans));  %
    spanLocs = bwlabel(~allInSpans);   %identify contiguous ones
    spikeRanges = unique(spanLocs); spikeRanges =spikeRanges(spikeRanges>0);
    meanRate = zeros(numel(spikeRanges),1);
    for jj = 1:numel(spikeRanges)
        meanRate(jj) = mean(spcnt(spanLocs==jj));
    end
    [~, id] = max(meanRate);
    zcnt = zscore(spcnt);
    goodZ = zcnt > -4 & zcnt < 4;
    goodSpikes = find(spanLocs==id & goodZ);
    
    
%     goodSpikes = find(spanLocs==1);
    plot(bins(goodSpikes), spcnt(goodSpikes), 'r.');
    if getInput
        tmp = input('enter spike range as a vector');
    else
        tmp = [];
    end
    if isempty(tmp) || ~(numel(tmp)==2)
        spikes.goodRange(units(ii),:) = [bins(goodSpikes(1)) bins(goodSpikes(end))];
    else
        spikes.goodRange(units(ii),:) = tmp;
    end
    xlabel('time (seconds)')
    ylabel('spike rate')
    title(sprintf('unit: %d', units(ii)))
end
