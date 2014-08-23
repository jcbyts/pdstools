function arOutput = getOutput(analysisRecord, outputName)
% GET analysis record output(s)
% arOutput = ov_existOutput(analysisRecord, outputName)

import ov.*

if ~isa(analysisRecord, 'us.physion.ovation.domain.concrete.AnalysisRecord')
	fprintf('first argument must be an ovation analysis record object\r')
	help getOutput
	return
end

outputs = ovation.map2struct(analysisRecord.getOutputs);

if isempty(outputs)
	arOutput = [];
	return
end

outputNames = fieldnames(outputs);
nOutputs = numel(outputNames);
if nargin < 2 || isempty(outputName)
    fprintf('found %d outputs\r', nOutputs)
    for ii = 1:nOutputs
        fprintf('%d. %s\r', ii, outputNames{ii})
    end
    arOutput = outputs;
    return
end

if isnumeric(outputName) && outputName <= nOutputs
    arOutput = outputs(outputName);
    return
end

isOutput = false(nOutputs,1);
for ii = 1:nOutputs
	isOutput(ii) = strcmp(outputName, outputNames{ii});
end

if any(isOutput)
	arOutput = outputs.(outputNames{isOutput});
else
	arOutput = [];
end