function spikes = getSpikes(pl, verbose, useContinuous)
% Get spike times and waveforms from plstruct 
% spikes = getSpikes(plstruct, verbose, useContinuousOnly)
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
if useContinuous
    switch pl.Version
        case 107 % offline sorter
            firstContinuousChannel = 64;
        case 106
            firstContinuousChannel = 64;

        otherwise
    end
else
    firstContinuousChannel = -1;
end
spikeChannels    = [pl.SpikeChannels(:).Channel]';
nUnitsPerChannel = zeros(numel(spikeChannels), 1);
for ii = 1:numel(spikeChannels)
    nUnitsPerChannel(ii) = numel(unique(pl.SpikeChannels(ii).Units));
end

channelsWithUnits = find(spikeChannels>firstContinuousChannel & nUnitsPerChannel>0);
nChannels = numel(channelsWithUnits);
spikeCtr = 0; % changed from 1 to 0... I believe we want to start with 0 because we add the unit number
spikes = struct();
spikes.time     = [];
spikes.id       = [];
spikes.waveform = [];
spikes.channel  = [];
spikes.timeaxis = [];
spikes.ADfreq   = pl.ADFrequency;

for ii = 1:nChannels
    spikes.time = [spikes.time; double(pl.SpikeChannels(channelsWithUnits(ii)).Timestamps)/pl.ADFrequency];
%     spikes.id   = [spikes.id; ones(numel(pl.SpikeChannels(channelsWithUnits(ii)).Timestamps),1).*double(pl.SpikeChannels(channelsWithUnits(ii)).Units(:)+spikeCtr)];
    spikes.id   = [spikes.id; double(pl.SpikeChannels(channelsWithUnits(ii)).Units(:))+spikeCtr];
    % convert spike waveforms to milliVolts

    % check if spike comes from AD continuous channel or sig channel
    SIG = pl.SpikeChannels(channelsWithUnits(ii)).SIGName;
    if strfind(SIG, 'sig')
        waves = (pl.SpikeMaxMagnitudeMV*double(pl.SpikeChannels(channelsWithUnits(ii)).Waves)') / (2^(pl.BitsPerSpikeSample)/2 * pl.SpikePreAmpGain * pl.SpikeChannels(channelsWithUnits(ii)).Gain);
    elseif strfind(SIG, 'adc')
        chMatch = [pl.ContinuousChannels(:).Channel] == pl.SpikeChannels(channelsWithUnits(ii)).Channel;
        if isempty(chMatch) || ~any(chMatch)
            waves = (pl.SpikeMaxMagnitudeMV*double(pl.SpikeChannels(channelsWithUnits(ii)).Waves)') / (2^(pl.BitsPerSpikeSample)/2 * pl.SpikePreAmpGain * pl.SpikeChannels(channelsWithUnits(ii)).Gain);
        else
            waves = (pl.ContMaxMagnitudeMV*double(pl.SpikeChannels(channelsWithUnits(ii)).Waves)') / (2^(pl.BitsPerContSample)/2 * pl.ContinuousChannels(chMatch).PreAmpGain * pl.ContinuousChannels(chMatch).ADGain);
        end
    else
        error('channel type not recognized. Gains not set.')
    end
        
    spikes.waveform = [spikes.waveform; waves];
%     spikes.waveform = [spikes.waveform; double(pl.SpikeChannels(channelsWithUnits(ii)).Waves)'];
%     spikes.channel  = [spikes.channel repmat(pl.SpikeChannels(channelsWithUnits(ii)).Channel, 1, pl.SpikeChannels(channelsWithUnits(ii)).NUnits)];
    spikes.channel  = [spikes.channel repmat(pl.SpikeChannels(channelsWithUnits(ii)).Channel, 1, nUnitsPerChannel(channelsWithUnits(ii)))];
%     spikeCtr = spikeCtr + pl.SpikeChannels(channelsWithUnits(ii)).NUnits;
     spikeCtr = spikeCtr + nUnitsPerChannel(channelsWithUnits(ii));
end

if useContinuous
    spikes.channel = spikes.channel - firstContinuousChannel;
end
[spikes.time, ord] = sort(spikes.time);
spikes.id          = spikes.id(ord);
spikes.waveform    = spikes.waveform(ord,:);
spikes.timeaxis    = ((1:pl.NumPointsWave) - pl.NumPointsPreThr)/pl.ADFrequency;
% calculate waveform SNR
units  = unique(spikes.id);
if any(units==0)
    spikes.id = spikes.id + 1;
    units = units + 1;
end

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
    units = unique(spikes.id);
    spn = ceil(sqrt(nUnits));
    for ii = 1:nUnits
        subplot(spn, spn, ii)
        idx = find(spikes.id==units(ii));
        plot(spikes.timeaxis*1e3, spikes.waveform(idx(1:100:end),:)', 'Color', .5*[1 1 1]); hold on
        plot(spikes.timeaxis*1e3, mean(spikes.waveform(idx, :)), 'k', 'Linewidth', 2); axis tight
        xlabel('msec')
        ylabel('mV')
        title(sprintf('un: %d, ch: %d, snr: %02.2f', ii, spikes.channel(ii), spikes.snr(ii)))
        
    end
    drawnow
end
