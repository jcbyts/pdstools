
%% testscript  for concatenateSpikes
import pdsa.*

load data/spikes1.mat
spikes1 = spikes;
load data/spikes2.mat
spikes2 = spikes;
clear spikes

[spikeStructs, spikesNew] = pdsa.concatenateSpikes(spikes1, spikes2)


%% testscript for readPLX

plxname = '/Volumes/LKCLAB/EPHYS/DATA/Projects/GaborDots/Pat/20130328/pat03282013gabordots1214-sort.plx';

pl = readPLXFileC(plxname);
% waveform time



%-----------------------------------------------------------------------------%
% make info struct
pl.WaveformTime = (-pl.NumPointsPreThr:(pl.NumPointsWave-pl.NumPointsPreThr))/pl.WaveformFreq;


%-----------------------------------------------------------------------------%
% make strobed struct
V = datevec(pl.Date);
assert(strcmp(pl.EventChannels(end).Name, 'Strobed'))
pl.Year = V(1);

strobeIndex = pl.EventChannels(end).Values == mod(pl.Year, 256);
strobed.times  = pl.EventChannels(end).Timestamps(strobeIndex);
strobed.values = reshape(pl.EventChannels(end).Values, 6, [])';

%-----------------------------------------------------------------------------%
% make events struct
nEvents = numel(pl.EventChannels)-1;
events.time = [];
events.id   = [];
events.name = cell(1);
eventCtr = 1;
for ii = 1:nEvents
	if pl.EventChannels(ii).Channel > 1 && pl.EventChannels(ii).Channel < 256
		events.time = [events.time; pl.EventChannels(ii).Timestamps];
		events.id   = [events.id; ones(numel(pl.EventChannels(ii).Timestamps),1)*eventCtr];
		events.name{eventCtr} = pl.EventChannels(ii).Name;
		eventCtr = eventCtr + 1;
	end
end

[events.time, ord] = sort(events.time);
events.id = events.id(ord);

%-----------------------------------------------------------------------------%
% make spikes struct
switch pl.Version
	case 107
		firstContinuousChannel = 64;
	otherwise

end
spikeChannels    = [pl.SpikeChannels(:).Channel]';
nUnitsPerChannel = [pl.SpikeChannels(:).NUnits]';

channelsWithUnits = find(spikeChannels>firstContinuousChannel & nUnitsPerChannel);
nChannels = numel(channelsWithUnits);
spikeCtr = 0;
spikes = struct();
spikes.time 	= [];
spikes.id   	= [];
spikes.waveform = [];
spikes.channel  = [];

for ii = 1:nChannels
	spikes.time = [spikes.time; double(pl.SpikeChannels(channelsWithUnits(ii)).Timestamps)/pl.ADFrequency];
	spikes.id   = [spikes.id; ones(numel(pl.SpikeChannels(channelsWithUnits(ii)).Timestamps),1).*double(pl.SpikeChannels(channelsWithUnits(ii)).Units(:)+spikeCtr)];
	spikes.waveform = [spikes.waveform; pl.SpikeChannels(channelsWithUnits(ii)).Waves'];
	spikes.channel  = [spikes.channel repmat(pl.SpikeChannels(channelsWithUnits(ii)).Channel, 1, pl.SpikeChannels(channelsWithUnits(ii)).NUnits)];
	spikeCtr = spikeCtr + 1;
end

[spikes.time, ord] = sort(spikes.time);
spikes.id 		   = spikes.id(ord);
spikes.waveform    = spikes.waveform(ord,:);

% calculate waveform SNR
nUnits = numel(unique(spikes.id));
spikes.snr = zeros(1,nUnits);
for ii = 1:nUnits
	spikes.id == 
end


