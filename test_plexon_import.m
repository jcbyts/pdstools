%% demo triggered average


% must be ABSOLUTE path to filename
plxname = '/Users/jacobyates/Dropbox/MatlabCode/Projects/2014_gaborDots/nancyChamberMappingMThyperflow/nancy20141023MT_L5P5_d10841.plx';

pl = readPLXFileC(plxname);

% get the sampling rate of each channel
samplingRates = [pl.ContinuousChannels(:).ADFrequency];

%-------------------------------------------------------------------------%
% there are *almost* always two sampling rates in a plexon file: a "fast"
% and a "slow" rate. "Fast" rate is the rate at which we sample the
% continuous high-pass filtered spike data. The "slow" rate is the sampling
% rate for LFP, Eyeposition and other analog data (emg?).

% find the slow rate, because that's where the LFP lives
slowRate = min(unique(samplingRates));

hasContinuousData = (pl.ContSampleCounts ~=0); % only take channels that have valid samples

continuousChannels = [pl.ContinuousChannels(samplingRates==slowRate & hasContinuousData).Channel]; % these are good LFP channels

nContinuousChannels = numel(continuousChannels);

%% Get Analog Data
tic
pl1 = readPLXFileC(plxname, 'continuous', continuousChannels);
toc

tic
[lfp_info, lfp_data] = plx.getAnalog(pl1, continuousChannels+1);
toc

pl2 = readPLXFileC(plname, 'events');


%% Get events
[events, strobed] = plx.getEvents(pl2);

%% Get Spikes
pls = readPLXFileC(plxname, 'spikes', 'waves');

spikes = plx.getSpikes(pls);

%% Trigger LFP on Events
ev = 1;
st = plx.convertTimeToSamples(events.time(events.id==ev), lfp_info.adfreq, lfp_info.timestamps, lfp_info.fragsamples);

[ta, ta_SD, bc] = pdsa.eventTriggeredAverage(lfp_data, st, [-500 500]);

plot(ta)
title(events.name{ev})

%% Trigger LFP on Spikes
un = 2;
st = plx.convertTimeToSamples(spikes.time(spikes.id==un), lfp_info.adfreq, lfp_info.timestamps, lfp_info.fragsamples);

[ta, ta_SD, bc] = pdsa.eventTriggeredAverage(lfp_data, st, [-500 500]);

figure(3); clf
plot(ta)
title(sprintf('units %d', un))

%% Trigger Spikes on Events
ev = 1;
un = 24;

[m,s,bc,v] = pdsa.eventPsth(spikes.time(spikes.id==un), events.time(events.id==ev), [-1 .1], 1e-3, ones(20,1)/20);

figure(4); clf
plot(bc, m)


%% plot Raster 

figure(4); clf
pdsa.plotRaster(spikes.time(spikes.id==un), events.time(events.id==ev), [-.1 1], 1e-3, [], 1)