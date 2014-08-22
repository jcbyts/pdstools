
%% testscript  for concatenateSpikes
import pdsa.*

load data/spikes1.mat
spikes1 = spikes;
load data/spikes2.mat
spikes2 = spikes;
clear spikes

[spikeStructs, spikesNew] = pdsa.concatenateSpikes(spikes1, spikes2);


%% testscript for readPLX

% plxname = '/Volumes/LKCLAB/EPHYS/DATA/Projects/GaborDots/Pat/20130328/pat03282013gabordots1214-sort.plx';
plxname = '/Users/jacobyates/Desktop/tmpPLXfiles/pat05012013MTmap1340-02.plx';


pl   = readPLXFileC(plxname, 'all');
%%
%-----------------------------------------------------------------------------%
% make info struct
info = readPLXFileC(plxname); % read header only
% waveform time
info.WaveformTime = (-pl.NumPointsPreThr:(pl.NumPointsWave-pl.NumPointsPreThr))/pl.WaveformFreq;


%-----------------------------------------------------------------------------%
% make strobed struct
[events, strobed] = plx.getEvents(pl);

%-----------------------------------------------------------------------------%
% make spikes struct
spikes = plx.getSpikes(pl);


%--------------------------------------------------------------------------------------------%
%% make analog struct:

if numel(pl.ContSampleCounts)~=numel(pl.ContinuousChannels)
    analogChannels = 1:numel(pl.ContinuousChannels);
else
    analogChannels = find(pl.ContSampleCounts>0);
end
analogSamplingRates  = [pl.ContinuousChannels(analogChannels).ADFrequency]';
analogRates = unique(analogSamplingRates); % to sort into lfp and analog
% if only one rate exists, only make the lfp struct
nRates = numel(analogRates);
if nRates > 0
	lfpChannels = analogChannels(analogSamplingRates == analogRates(1));
	[lfp_info, lfp_data] = plx.getAnalog(pl, lfpChannels);
	if nRates > 1
		contChannels = analogChannels(analogSamplingRates == analogRates(2));
		[an_info, an_data] = plx.getAnalog(pl, contChannels);
	end
end