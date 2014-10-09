function epoch = getOrInsertEpoch(experiment, protocol, plxfile, dv, Timezone)
% GET or ADD epoch from an experiment object
% epoch = getOrInsertEpoch(experiment, protocol_name, varargin)
% ov.getOrInsertEpoch can add or retrieve an epoch. Returns a boolean if the 
% specified epoch does not exist and the required variables were not proveided
% to create it.
% epoch = ov.getOrInsertEpoch(experiment, protocol_name)
% 	experiment  [java obj] 
% 	protocol_name [string] - name of the ptrotocl
% returns epoch object if the epoch exists. returns false if it doesn't.
% epoch = ov.getOrInsertEpoch(experiment, protocol_name, plxfile, dv)
% 	experiment  [java obj] 
% 	protocol_name [string] - name of the ptrotocl
% 		  plxfile [string] - full path to plexon file
% 		       dv [struct] - matlab pldaps struct
% 	returns epoch if it exists. Adds it if it doesn't. Name it takes
%  [protocol_name]. startDate and endDate are calculated from the plxfile.
%  parameters and device parameters are added with dv.pa and dv.disp
% 
% dependencies:
% readPLXFileC

import ov.*

if ~isa(experiment, 'us.physion.ovation.domain.concrete.Experiment')
	fprintf(['\n\n*****************************************************************\n' ...
		'*****************************************************************\n\n\n']);
	fprintf('You need to pass in an ovation Experiment object as the first argument.\nYou passed in a %s object\n\n', class(experiment))
	if isa(experiment, 'us.physion.ovation.api.StdDataContext')
		experiment = ov.getExperiment(experiment, [], []);
	elseif isa(experiment, 'us.physion.ovation.domain.concrete.Project')
		experiment = ov.getExperiment(experiment, []);
	end
	return
end

epochs = ovation.asarray(experiment.getEpochs);
nEpochs = numel(epochs);

if nargin < 2
	protocol_name = '';
else
	if ischar(protocol)
		protocol_name = protocol;
	elseif isnumeric(protocol)
		if protocol <=nEpochs
            ep = epochs(protocol);
            prot = ep.getProtocol;
            if isempty(prot)
                epoch = [];
                return
            end
			protocol_name = char(prot.getName);
		else
			epoch = epochs;
			return
		end
	else
		protocol_name = char(protocol.getName);
	end
end

if ~exist('Timezone', 'var')
    Timezone = 'America/Chicago';
end
% 20140614 jly 	wrote it





epochnames  = cell(nEpochs,1);
iscondition = false(nEpochs,1);
for ii = 1:nEpochs
	p = epochs(ii).getProtocol;
	if ~isempty(p)
		epochnames{ii} = char(p.getName);
		iscondition(ii) = strcmp(epochnames{ii}, protocol_name);
	end
end

if nargin == 1
	fprintf('found %d epochs\n', nEpochs)
	for ii = 1:nEpochs
		fprintf('%s\n', epochnames{ii})
	end
	epoch = epochs;
	return

end



if any(iscondition)
	epoch = epochs(iscondition);
    return
elseif nargin == 2
	epoch = [];
    return
end

%-----------------------------------------------------------------------%
% adding epoch group
pl = readPLXFileC(plxfile);
[year, month, day, hour, minute, seconds] = datevec(pl.Date);

startDate = ovation.datetime(year, month, day, hour, minutes, seconds,0, Timezone);
endTime.full = hour*60*60+minutes*60+seconds+(pl.LastTimestamp/pl.WaveformFreq);
endTime.hour = floor(endTime.full/60/60);
endTime.minute = floor((endTime.full-endTime.hour*60*60)/60);
endTime.second = floor(endTime.full-(endTime.hour*60*60+endTime.minute*60));
endDate   = ovation.datetime(year, month, day, endTime.hour,endTime.minute, endTime.second, 0, Timezone);
% guess end time
% [~, adfreq] = plx_adchan_freqs(plxfile);
% [~, adc]    = plx_adchan_samplecounts(plxfile);
if iscell(dv)
    tmpdv = dv{40};
else
    tmpdv = dv;
end
epoch = experiment.insertEpoch(startDate, endDate, protocol, ovation.struct2map(tmpdv.pa), ovation.struct2map(tmpdv.disp));
return
