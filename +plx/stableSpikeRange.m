function spikes = stableSpikeRange(spikes)
% find range with stable recordings

units = unique(spikes.id);
nUnits = numel(units);

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
   figure(1); clf; hold all
    plot(bins, spcnt, 'o', bins, rhat, 'r')
%     plot(bins, zscore(spcnt), 'o')

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
    goodSpikes = find(spanLocs==id);
%     goodSpikes = find(spanLocs==1);
    plot(bins(goodSpikes), spcnt(goodSpikes), 'r.');
    spikes.goodRange(units(ii),:) = [bins(goodSpikes(1)) bins(goodSpikes(end))];

    pause
end