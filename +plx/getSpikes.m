function spikes = plx_getSpikes(plxname, useContinuous, info, verbose)
% spikes = plx_getSpikes(plxname, continuousOnly, info, verbose)
% get spike times and waveforms from plxname.plx
% inputs:
%   plxname     [string] - full/path/to/*.plx
%   continuousOnly[bool] - if 1, use only spikes that came from a channel
%                          with analog sample. This is required for binary
%                          pursuit.
%                          if 0, get all spikes and waveforms from any
%                          channel with spikes. If you sorted both spike
%                          channels and 'spike+continuous' channels there
%                          will be duplicates of units. This setting should
%                          be used if you didn't sort continuous channels.
%   info        [struct] - meta info about the experiment
%       .year - required to parse strobe values

if nargin < 4
    verbose = 1;
    if nargin < 3
        info = [];
        if nargin < 2
            useContinuous = 0;
            if nargin < 1
                [plxfile, plxpath] = uigetfile('*.plx');
                plxname = fullfile(plxpath, plxfile);
            end
        end
    end
end

if isempty(info) || ~isfield(info, 'waveform_time')
    info = plx_getInfo(plxname, info);
end

spikes = [];

%-------------------------------------------------------------------------%
% find channels that have spike events on them
[~, wfcounts] = plx_info(plxname,1);

[spike_unit,spike_channels]  = find(wfcounts);
[~,dspchans]                 = plx_chanmap(plxname);
% this includes all channels that are classified as 'spike' channels in
% plexon's offline sorter.
spike_channels               = dspchans(spike_channels-1); % channels that could contain spikes


%-------------------------------------------------------------------------%
% find all channels that have analog samples
% get channel map
[~,sdc] = plx_ad_chanmap(plxname);
% [~, chan_names] = plx_adchan_names(plxname); % maybe use this to remove eye channels
[~,freqs]   = plx_adchan_freqs(plxname);
Fs = unique(freqs); % sampling rates in the file

continuous_channels = sdc(freqs==max(Fs));
% correct spike channels for any units that were sorted on continuous
% channels
first_continuous = continuous_channels(1);

nSpkCh = length(spike_channels);
if useContinuous
    if verbose
        fprintf('finding %d units on %d channels\r', sum(spike_channels>first_continuous), numel(unique(spike_channels(spike_channels>first_continuous))))
    end
else
    if verbose
        fprintf('finding %d units on %d channels\r', nSpkCh, numel(unique(spike_channels)))
    end
end
spike_waveforms  = []; % or preallocate? ... sum(wfcounts(:))
spike_times = [];
spike_id         = [];
snr = double(nSpkCh);
unctr =1;
for ii = 1:nSpkCh
    if useContinuous && spike_channels(ii) < first_continuous
        spike_channels(ii) = 0;
    else
        if verbose
            fprintf('loading channel: %0.0i, unit: %0.0i \r', spike_channels(ii), spike_unit(ii))
        end
        [n,~, wave_ts, spike_waves] = plx_waves(plxname, spike_channels(ii), spike_unit(ii)-1);
        spike_waveforms  = [spike_waveforms; spike_waves]; %#ok<AGROW>
        spike_times = [spike_times; wave_ts]; %#ok<AGROW>
        spike_id         = [spike_id; ones(n,1)*unctr]; %#ok<AGROW>
        % spike waveform signal-to-noise... measure of single unitness
        % for equations see Kelly et al 2007 (Movshon is final author)
        mw = mean(spike_waves);
        peak   = max(mw);
        trough = min(mw);
        amp = peak - trough;
        mwbar = repmat(mw, n, 1);
        r = spike_waves - mwbar;
        sigma = std(r(:));
        snr(unctr) = amp/(2*sigma);
        unctr = unctr+1;
    end
end
spike_channels(spike_channels>first_continuous) = spike_channels(spike_channels>first_continuous);
spike_channels(spike_channels==0) = [];
if verbose
    disp('Done!')
end
nUnits = numel(snr);
%-------------------------------------------------------------------------%
% plot waveforms
if verbose
    figure(1); clf
    spn = ceil(sqrt(nUnits));
    for ii = 1:nUnits
        subplot(spn, spn, ii)
        idx = find(spike_id==ii);
        plot(info.waveform_time*1e3, spike_waveforms(idx(1:100:end),:), 'Color', .5*[1 1 1]); hold on
        plot(info.waveform_time*1e3, mean(spike_waveforms(idx, :)), 'k', 'Linewidth', 2); axis tight
        xlabel('time')
        ylabel('micro volts?')
        title(sprintf('un: %d, ch: %d, snr: %02.2f', ii, spike_channels(ii), snr(ii)))
        
    end
    drawnow
end

[sptime, sortid] = sort(spike_times);
spikes.time     = sptime;
spikes.id       = spike_id(sortid);
spikes.waveform = spike_waveforms(sortid,:);
spikes.channel  = spike_channels;
spikes.snr      = snr;
spikes.first_continuous_channel = first_continuous;
spikes.continuous_only = useContinuous;