function plotSpikeWaveforms(spikes, SNRthresh, Chans)
% plot spike waveforms by channel
% plotSpikeWaveforms(spikes, SNRthresh, Chans)
% Inputs:
%   spikes [struct]
%       .waveform
%       .id
%       .snr
%       .channel
%   SNRthresh [1 x 1] - threshold for considering single units (default=3)
%   Chans     [1 x m] - vector of channels to consider

% plotSpikes
if ~exist('SNRthresh', 'var')
    SNRthresh = 3;
end

if ~exist('Chans', 'var')
    Chans = 1:max(spikes.channel);
end

% nUnits = numel(spikes.snr);
% nUnits = sum(sum(bsxfun(@eq, spikes.channel, Chans'),1)>0);
nUnits = sum(spikes.snr(sum(bsxfun(@eq, spikes.channel, Chans'),1)>0)>SNRthresh);
cmap = hsv(nUnits);

hold all
nSamples = size(spikes.waveform,2);

k = 1;
channels = unique(spikes.channel);
channels = intersect(channels, Chans);
for ii = 1:numel(channels)
    unitsOnChannel = find(spikes.channel == channels(ii));
    [~, ord] = sort(spikes.snr(unitsOnChannel));
    
    
    for id = unitsOnChannel(ord)
        if spikes.snr(id) < SNRthresh
            co = .5*[1 1 1];
        else
            co = cmap(k,:);
        end
            Unit = find(spikes.id==id);
            n = numel(Unit);
            s = ceil(n/100);
            wave = spikes.waveform(Unit(1:s:end),:)./abs(min(min(spikes.waveform(Unit(1:s:end),:))));
            plot(1:nSamples, wave-2*channels(ii), 'Color', co);
            k = k+1;
    end
end

set(gca, 'YTick', sort(-2*channels), 'YTickLabel', abs(sort(-channels)))
axis tight
ylabel('Channel')
