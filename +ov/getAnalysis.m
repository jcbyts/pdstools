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

% assert(isa(epoch, 'us.physion.ovation.domain.concrete.Epoch'), 'looking for analysis records for a specific epoch')

anRecs = ovation.asarray(epoch.getAnalysisRecords());
nAnRecs = numel(anRecs);

if nAnRecs == 0
	analysisRecord = [];
	fprintf('there are no analysis records for this epoch\n')
	return
end

if nargin < 2 || isempty(analysisRecordName)
	fprintf('found %d analysis records.\n', nAnRecs)
	for ii = 1:nAnRecs
		fprintf('%d. %s\n', ii, char(anRecs(ii).getName()))
	end
	analysisRecord = anRecs;
	return
end

if isnumeric(analysisRecordName) && analysisRecordName<=nAnRecs
	analysisRecord = anRecs(analysisRecordName);
	fprintf('returning analysis record: %s\n', char(analysisRecord.getName))
	return
end


if nAnRecs>0
	isAnal = false(nAnRecs,1);
	if isa(analysisRecordName, 'us.physion.ovation.domain.concrete.AnalysisRecord')
		for ii = 1:nAnRecs
			isAnal(ii) = analysisRecordName == anRecs(ii);
		end
	elseif ischar(analysisRecordName)
		for ii = 1:nAnRecs
			isAnal(ii) = strcmp(analysisRecordName, anRecs(ii).getName());
		end
	end
else
	isAnal = false;
end

if any(isAnal)
	analysisRecord = anRecs(isAnal);
else
	analysisRecord = [];
end
