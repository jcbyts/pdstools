function importPLXEvents(epoch, pl)
% import plexon events and strobed into ovation analysis record
% importPLXEvents(epoch, pl)

import ov.*

assert(isa(epoch, 'us.physion.ovation.domain.concrete.Epoch'), 'first argument must be an epoch')

rawPLX 	= getAnalysis(epoch, 'Raw PLX Data');


if nargin < 2 || isempty(pl)
	[plxfile, plxpath] = uigetfile('*.plx');
	pl = fullfile(plxpath, plxfile);
end

if isstr(pl)
	if exist(pl, 'file')
		fprintf('reading in plexon file')
		pl   = readPLXFileC(pl, 'all');
	else
		fprintf('fail\r')
		return
	end
end

[events, strobed] = plx.getEvents(pl);

spikes = plx.getSpikes(pl);


infoFields = fieldnames(pl);
info = struct();
for ii = 1:21
 	info.(infoFields{ii}) = pl.(infoFields{ii});
end

info.WaveformTime = (-pl.NumPointsPreThr:(pl.NumPointsWave-pl.NumPointsPreThr))/pl.WaveformFreq;

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



if isempty(rawPLX)
	epoch.addAnalysisRecord('Raw PLX Data', epoch.getMeasurements, epoch.getProtocol, ovation.struct2map(info));
	rawPLX 	= getAnalysis(epoch, 'Raw PLX Data');
end

if isempty(ov.getOutput(rawPLX, 'continuous'))
	mkdir('tmp')
	fname = [char(epoch.getExperiment.getPurpose) '_continuous.mat'];
	fname = fullfile(pwd, 'tmp', fname);
	save(fname, 'an_info', 'an_data')
	addOutput(rawPLX, fname, 'continuous')
	delete(fname)
	rmdir('tmp')
end

if isempty(ov.getOutput(rawPLX, 'lfp'))
	mkdir('tmp')
	fname = [char(epoch.getExperiment.getPurpose) '_lfp.mat'];
	fname = fullfile(pwd, 'tmp', fname);
	save(fname, 'lfp_info', 'lfp_data')
	addOutput(rawPLX, fname, 'lfp')
	delete(fname)
	rmdir('tmp')
end

if isempty(ov.getOutput(rawPLX, 'events'))
	mkdir('tmp')
	fname = [char(epoch.getExperiment.getPurpose) '_events.mat'];
	fname = fullfile(pwd, 'tmp', fname);
	save(fname, 'events', 'strobed', 'info')
	addOutput(rawPLX, fname, 'events')
	delete(fname)
	rmdir('tmp')
end

if isempty(ov.getOutput(rawPLX, 'spikes'))
	mkdir('tmp')
	fname = [char(epoch.getExperiment.getPurpose) '_spikes.mat'];
	fname = fullfile(pwd, 'tmp', fname);
	save(fname, 'spikes')
	addOutput(rawPLX, fname, 'spikes')
	delete(fname)
	rmdir('tmp')
end
