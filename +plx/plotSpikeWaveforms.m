function plotSpikeWaveforms(spikes, SNRthresh, Chans, MaxChannel)
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
%   MaxChannel[1 x 1] - maximum channel on rig (loops channels above that)

% plotSpikes
if ~exist('SNRthresh', 'var')
    SNRthresh = 3;
end

if ~exist('MaxChannel', 'var')
    MaxChannel = max(spikes.channel);
end

if ~exist('Chans', 'var') || isempty(Chans)
    Chans = 1:max(spikes.channel);
end

% nUnits = numel(spikes.snr);
% nUnits = sum(sum(bsxfun(@eq, spikes.channel, Chans'),1)>0);
nUnits = sum(spikes.snr(sum(bsxfun(@eq, spikes.channel, Chans'),1)>0)>SNRthresh);
cmap = hsv(nUnits);

nSamples = size(spikes.waveform,2);
% set x-axis to time if it exists
if isfield(spikes, 'timeaxis')
    xax = spikes.timeaxis*1e3; % convert to milliseconds
    units = 'msec';
else
    xax = 1:nSamples;
    units = 'samples';
end

hold all


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
            k = k+1;
        end
            Unit = find(spikes.id==id);
            n = numel(Unit);
            s = ceil(n/100);
            wave = spikes.waveform(Unit(1:s:end),:)./abs(min(min(spikes.waveform(Unit(1:s:end),:))));
            plot(xax, wave-2*mod(channels(ii), MaxChannel), 'Color', co);
    end
end

set(gca, 'YTick', sort(-2*channels), 'YTickLabel', abs(sort(-channels)))
axis tight
ylabel('Channel')
xlabel(units)
