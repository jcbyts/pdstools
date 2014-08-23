function analysisRecord = getAnalysis(epoch, analysisRecordName)
% GET analysis record(s) from epoch
% analysisRecord = getAnalysis(epoch, analysisRecordName)
% Check if analysis record exists
% INPUTS
% 	epoch [object] - ovation epoch object
% 	analysisRecordName - multiple ways to call
% 					   [string] - checks if string exists as an analysis. 
% 								  Returns empty if it doesn't.
% 					   [object]	- checks if an object exists as an anlaysis 
% 								  record. Returns empty if it doesn't.
% 					    [empty] - prints list of existing analysis records
%  					  [numeric] - returns the [numeric]th  analysis record.
% 								Use ov_existAnalysis(epoch, []) to see list of 
% 								available analysis records

import ov.*

assert(isa(epoch, 'us.physion.ovation.domain.concrete.Epoch'), 'looking for analysis records for a specific epoch')

anals = ovation.asarray(epoch.getAnalysisRecords());
nAnals = numel(anals);

if nAnals == 0
	analysisRecord = [];
	fprintf('there are no analysis records for this epoch\r')
	return
end

if nargin < 2 || isempty(analysisRecordName)
	fprintf('found %d analysis records.\r', nAnals)
	for ii = 1:nAnals
		fprintf('%d. %s\r', ii, char(anals(ii).getName()))
	end
	analysisRecord = anals;
	return
end

if isnumeric(analysisRecordName) && analysisRecordName<=nAnals
	analysisRecord = anals(analysisRecordName);
	fprintf('returning analysis record: %s\r', char(analysisRecord.getName))
	return
end


if nAnals>0
	isAnal = false(nAnals,1);
	if isa(analysisRecordName, 'us.physion.ovation.domain.concrete.AnalysisRecord')
		for ii = 1:nAnals
			isAnal(ii) = analysisRecordName == anals(ii);
		end
	elseif ischar(analysisRecordName)
		for ii = 1:nAnals
			isAnal(ii) = strcmp(analysisRecordName, anals(ii).getName());
		end
	end
else
	isAnal = false;
end

if any(isAnal)
	analysisRecord = anals(isAnal);
else
	analysisRecord = [];
end