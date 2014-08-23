function spikes = getSpikes(pl, verbose, useContinuous)
% spikes = getSpikes(plstruct, verbose, useContinuousOnly)
% get spike times and waveforms from plstruct
% inputs:
%   plstruct    [struct] - output of readPLXFileC(plxname, 'all')
%   verbose     [0 or 1] - plot stuff if 1
%   continuousOnly[bool] - if 1, use only spikes that came from a channel
%                          with analog sample. This is required for binary
%                          pursuit.
%                          if 0, get all spikes and waveforms from any
%                          channel with spikes. If you sorted both spike
%                          channels and 'spike+continuous' channels there
%                          will be duplicates of units. This setting should
%                          be used if you didn't sort continuous channels.
% outputs       
%   spikes  [struct]

spikes = [];

if nargin < 3
    verbose = 1;
    if nargin < 2
        useContinuous = 0;
        if nargin < 1
            help plx.getSpikes
            return
        end
    end
end


%-----------------------------------------------------------------------------%
%% make spikes struct
% this will be different for each rig. I want to find a nice way to handle it, but we
% don't currently have one. It might hav eto be a manual option.
switch pl.Version
    case 107 % offline sorter
        firstContinuousChannel = 64;
    case 106
        firstContinuousChannel = 64;

    otherwise

end
spikeChannels    = [pl.SpikeChannels(:).Channel]';
nUnitsPerChannel = [pl.SpikeChannels(:).NUnits]';

channelsWithUnits = find(spikeChannels>firstContinuousChannel & nUnitsPerChannel);
nChannels = numel(channelsWithUnits);
spikeCtr = 0;
spikes = struct();
spikes.time     = [];
spikes.id       = [];
spikes.waveform = [];
spikes.channel  = [];

for ii = 1:nChannels
    spikes.time = [spikes.time; double(pl.SpikeChannels(channelsWithUnits(ii)).Timestamps)/pl.ADFrequency];
    spikes.id   = [spikes.id; ones(numel(pl.SpikeChannels(channelsWithUnits(ii)).Timestamps),1).*double(pl.SpikeChannels(channelsWithUnits(ii)).Units(:)+spikeCtr)];
    spikes.waveform = [spikes.waveform; double(pl.SpikeChannels(channelsWithUnits(ii)).Waves)'];
    spikes.channel  = [spikes.channel repmat(pl.SpikeChannels(channelsWithUnits(ii)).Channel, 1, pl.SpikeChannels(channelsWithUnits(ii)).NUnits)];
    spikeCtr = spikeCtr + 1;
end

[spikes.time, ord] = sort(spikes.time);
spikes.id          = spikes.id(ord);
spikes.waveform    = spikes.waveform(ord,:);

% calculate waveform SNR
units  = unique(spikes.id);
nUnits = numel(units);
spikes.snr = zeros(1,nUnits);
for ii = 1:nUnits
    idx = spikes.id == units(ii);
    mw = mean(spikes.waveform(idx,:));
    peak   = max(mw);
    trough = min(mw);
    amp    = peak - trough;
    r  = spikes.waveform(idx,:) - repmat(mw, sum(idx),1);
    sigma  = std(r(:));
    spikes.snr(ii) = amp/(2*sigma);
end

%-------------------------------------------------------------------------%
% plot waveforms
if verbose
    figure(1); clf
    spn = ceil(sqrt(nUnits));
    for ii = 1:nUnits
        subplot(spn, spn, ii)
        idx = find(spikes.id==ii);
        plot(spikes.waveform(idx(1:100:end),:)', 'Color', .5*[1 1 1]); hold on
        plot(mean(spikes.waveform(idx, :)), 'k', 'Linewidth', 2); axis tight
        xlabel('samples')
        ylabel('micro volts?')
        title(sprintf('un: %d, ch: %d, snr: %02.2f', ii, spikes.channel(ii), spikes.snr(ii)))
        
    end
    drawnow
end