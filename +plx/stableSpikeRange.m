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
        binSize = 30; % seconds
        bins = 0:binSize:max(stimes);
        spcnt = histc(stimes, bins)/binSize;
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
    
    binSize = 60; % seconds
    timeBins = 0:binSize:max(stimes);
    binnedSpikeRate = histc(stimes, timeBins)/binSize;
    
    % eliminate times when the unit did not fire any spikes in the window
    goodIndex=binnedSpikeRate~=0;
    
%     bins=timeBins(goodIndex);
%     spcnt=binnedSpikeRate(goodIndex);
    
    % fit line
    x = [timeBins(:) ones(numel(timeBins),1)];
    wls = (x(goodIndex,:)'*x(goodIndex,:))\(x(goodIndex,:)'*binnedSpikeRate(goodIndex));
%     lm=fitlm(timeBins(goodIndex), binnedSpikeRate(goodIndex));
%     
%     x = [bins(:) ones(numel(bins),1)];
%     wls = (x'*x)\(x'*spcnt);
    rhat = x*wls;
    residuals = binnedSpikeRate - rhat; 
%     resids = spcnt - rhat; 
    h(ii) = figure(ii); clf; hold all
	plot(timeBins, binnedSpikeRate, 'o', timeBins, rhat, 'r')

    % operate only on the good indexed spikes
    residualsGood = residuals(goodIndex);
    spikeRateGood   = binnedSpikeRate(goodIndex);
    % find regions where the line is a good fit
    
    belowThreshold = abs(residualsGood/std(residualsGood))<3;
    spanLocs = bwlabel(belowThreshold);   %identify contiguous ones
    spanLength = regionprops(spanLocs, 'area');  %#ok<MRPBW> %length of each span
    spanLength = [spanLength.Area];
    
    [~, id]=max(spanLength); % largest region of stability
    splist=find(goodIndex);
    goodSpikes = splist(spanLocs==id);
    
%     goodSpans = find(spanLength>=5);   %get on ly spans of 5+ points
%     allInSpans = (ismember(spanLocs, goodSpans));
%     spanLocs = bwlabel(allInSpans);   %identify contiguous ones
%     spikeRanges = unique(spanLocs); spikeRanges=spikeRanges(spikeRanges>0);
%     meanRate = zeros(numel(spikeRanges),1);
%     for jj = 1:numel(spikeRanges)
%         meanRate(jj) = mean(spikeRateGood(spanLocs==jj));
%     end
%     [~, id] = max(meanRate);
%     splist=find(goodIndex);
%     goodSpikes = splist(spanLocs==id);
    
%     % extra check for large fluctuations
%     zcnt = zscore(spikeRateGood);
%     goodZ = zcnt > -4 & zcnt < 4;
%     
%     goodSpikes = splist(find(spanLocs==id & goodZ));
    
    
%     belowThreshold = abs(residuals/std(residuals))<3;
%     spanLocs = bwlabel(belowThreshold);   %identify contiguous ones
%     spanLength = regionprops(spanLocs, 'area');  %#ok<MRPBW> %length of each span
%     spanLength = [ spanLength.Area];
%     goodSpans = find(spanLength>=5);   %get on ly spans of 5+ points
%     allInSpans = (ismember(spanLocs, goodSpans));  %
%     spanLocs = bwlabel(allInSpans);   %identify contiguous ones
%     spikeRanges = unique(spanLocs); spikeRanges =spikeRanges(spikeRanges>0);
%     meanRate = zeros(numel(spikeRanges),1);
%     for jj = 1:numel(spikeRanges)
%         meanRate(jj) = mean(spcnt(spanLocs==jj));
%     end
%     [~, id] = max(meanRate);
%     zcnt = zscore(spcnt);
%     goodZ = zcnt > -4 & zcnt < 4;
%     splist=1:numel(timeBins);
%     goodSpikes = splist(find(spanLocs==id & goodZ));
    
    
%     goodSpikes = find(spanLocs==1);
%     plot(bins(goodSpikes), spcnt(goodSpikes), 'r.');
    plot(timeBins(goodSpikes), binnedSpikeRate(goodSpikes), 'r+');
    if getInput
        tmp = input('enter spike range as a vector');
    else
        tmp = [];
    end
    if isempty(tmp) || ~(numel(tmp)==2)
        spikes.goodRange(units(ii),:) = [timeBins(goodSpikes(1)) timeBins(goodSpikes(end))];
    else
        spikes.goodRange(units(ii),:) = tmp;
    end
    xlabel('time (seconds)')
    ylabel('spike rate')
    title(sprintf('unit: %d', units(ii)))
end
