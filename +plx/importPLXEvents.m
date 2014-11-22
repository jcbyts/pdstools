function importPLXEvents(epoch, pl)
% import plexon events and strobed into ovation analysis record
% importPLXEvents(epoch, pl)
% This needs serious documentation work
import ov.*
import plx.*
import pdsa.*

assert(isa(epoch, 'us.physion.ovation.domain.concrete.Epoch'), 'first argument must be an epoch')

rawPLX 	= getAnalysis(epoch, 'Raw PLX Data');


if nargin < 2 || isempty(pl)
	meas = ov.getOrInsertMeasurement(epoch, 'plx');
	nMeas = numel(meas);
	if nMeas == 1
		pl = ovation.datapath(meas);
	elseif nMeas > 1
		for ii = 1:nMeas
			fprintf('%d) %s\n', ii, char(meas(ii).getName))
		end
		mi = input('which plx file did you want to load?'  );
		pl = ovation.datapath(meas(mi));
	else
		[plxfile, plxpath] = uigetfile('*.plx');
		pl = fullfile(plxpath, plxfile);
	end
end

if isstr(pl)
	if exist(pl, 'file')
		fprintf('reading in plexon file\n')
		pl   = readPLXFileC(pl, 'all');
	else
		fprintf('fail\n')
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

if numel(rawPLX) > 1
    rawPLX = rawPLX(1);
end

if isempty(ov.getOutput(rawPLX, 'continuous')) && exist('an_info', 'var')
	mkdir('tmp')
	fname = [char(epoch.getExperiment.getPurpose) '_continuous.mat'];
	fname = fullfile(pwd, 'tmp', fname);
	save(fname, 'an_info', 'an_data', '-v7.3')
    fprintf('adding continuous data\n')
	addOutput(rawPLX, fname, 'continuous')
	delete(fname)
	rmdir('tmp')
    fprintf('Done\n')
end

if isempty(ov.getOutput(rawPLX, 'lfp')) && exist('lfp_info', 'var')
	mkdir('tmp')
	fname = [char(epoch.getExperiment.getPurpose) '_lfp.mat'];
	fname = fullfile(pwd, 'tmp', fname);
	save(fname, 'lfp_info', 'lfp_data', '-v7.3')
    fprintf('adding lfp data\n')
	addOutput(rawPLX, fname, 'lfp')
	delete(fname)
	rmdir('tmp')
    fprintf('Done\n')
end

if isempty(ov.getOutput(rawPLX, 'events'))
	mkdir('tmp')
	fname = [char(epoch.getExperiment.getPurpose) '_events.mat'];
	fname = fullfile(pwd, 'tmp', fname);
	save(fname, 'events', 'strobed', 'info')
    fprintf('adding event data\n')
	addOutput(rawPLX, fname, 'events')
	delete(fname)
	rmdir('tmp')
    fprintf('Done\n')
end

if isempty(ov.getOutput(rawPLX, 'spikes'))
	mkdir('tmp')
	fname = [char(epoch.getExperiment.getPurpose) '_spikes.mat'];
	fname = fullfile(pwd, 'tmp', fname);
    fprintf('adding spike data\n')
	save(fname, 'spikes')
	addOutput(rawPLX, fname, 'spikes')
	delete(fname)
	rmdir('tmp')
    fprintf('Done\n')
end
