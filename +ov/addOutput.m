function addOutput(analysisRecord, saveFileName, outputName) 
% ADD output file to an ovation analysis record object
% addOutput(analysisRecord, saveFileName, outputName)
% inputs:
% 	analysisRecord  [obj] - analysis record object
% 	saveFileName [string] - full/path/to/file
% 	outputName   [string] - what to call the output

% 20140531 jly	wrote it
% 20140822 jly	package version

assert(isa(analysisRecord, 'us.physion.ovation.domain.concrete.AnalysisRecord'), 'add output only works for analysis records')
assert(exist(saveFileName, 'file')>0, 'saveFileName must be a path to a file that exists')
assert(ischar(outputName), 'you can call your output whatever you want, but it has to be a string')

f = java.io.File(saveFileName);
if f.isAbsolute
		url = f.toURI().toURL();
		analysisRecord.addOutput(outputName, url, ovation.util.content_type(saveFileName));
else
	warning('filename has to be an absolute path!')
	return
end
