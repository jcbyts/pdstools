function measurement = getOrInsertMeasurement(epoch, tag, files)
% GET or ADD measurement from/to an epoch
% measurement = getOrInsertMeasurement(epoch, tag, files)
%
%	measurement = ov_getOrInsertMeasurement(epoch, tag);
% 		gets measurements that have the tag [tag] from 
% 		an ovation epoch. 
%   measurement = ov_getOrInsertMeasurement(epoch, tag, files)
% 		returns measurements that match tag and files. If the file doesn't
% 		exist, it gets added with the tag. If the tag doesn't exist, returns empty.
%
% Inputs:
% 	epoch [java obj] - ovation epoch object
% 		tag [string] - string for the tag
% files [cell array] - cellarray of strings with absolute file path
%
% Output:
% 	measurement - array of measurement objects

% 20140614 jly 	wrote it

if nargin < 3
	files = [];
	if nargin < 2
		tag = [];
	end
end

measurements = ovation.asarray(epoch.getMeasurements);
nMeasurements = numel(measurements);


if isempty(files)
if isempty(tag)
	measurement = measurements;
	fprintf('found %d measurements\n', nMeasurements)
	for ii = 1:nMeasurements
		fprintf('%d) %s\n', ii, char(measurements(ii).getName))
	end
	return
elseif isnumeric(tag) && tag <= nMeasurements
	measurement = measurements(tag);
	return
elseif ischar(tag)
	measBool = false(nMeasurements,1);
	for ii = 1:nMeasurements
        tags = ovation.asarray(measurements(ii).getAllTags);
        nTags = numel(tags);
        tagMatch = false(numel(tags),1);
        for jj = 1:nTags
            tagMatch(jj) = strcmp(char(tags(1)), tag);
        end
        measBool(ii) = any(tagMatch);        
	end		
	fprintf('found %d measurements that match %s\n', sum(measBool), tag)
	measurement = measurements(measBool);
	return

end

else
	body

	if ischar(files)
		files = {files};
	end
	nFiles = numel(files);
	fileInd = [];
	for ii = 1:nFiles
		[measurement, id] = ov_hasFile(measurements, files{ii});

		% if the file doesn't exist, add it with the tag
		if isempty(measurement)
			[pt, fn, ext] = fileparts(files{ii});
			f = java.io.File(files{ii});
			if f.isAbsolute
				url = f.toURI().toURL();
				measurement = epoch.insertMeasurement([fn ext], [],[], url, 'application/x-pldaps');
				if ischar(tag)
					measurement.addTag(tag);
				end
				% get measurements again and check for file, it should
				% be there now
				measurements = ovation.asarray(epoch.getMeasurements());
				[measurement, id] = ov_hasFile(measurements, files{ii});
			end
		end
		% iterate to get file indices
		fileInd = [fileInd; id];
	end
	% return only the measurements that match the file
	measurements = measurements(fileInd);
end

% if ~isempty(measurements)
% 	measurement = ov_getTag(measurements, tag);
% end

end


function [measurement, id] = ov_hasFile(measurements, file)
	if isempty(measurements)
		measurement = [];
		id = [];
		return
	end

	nMeasurements = numel(measurements);
	[pt, fn, ext] = fileparts(file);
	fileBool = false(nMeasurements,1);
	for jj = 1:nMeasurements
		fileBool(jj) = strcmp([fn ext], char(measurements(jj).getName));
	end
	if any(fileBool)
		measurement = measurements(fileBool);
		id = find(fileBool);
	else
		measurement = [];
		id = [];
	end
end


function measurement = ov_getTag(measurements, tag)
	if isempty(measurements)
		measurement = [];
		return
	end

	if isempty(tag)
		measurement = measurements;
		return
	end

	nMeasurements = numel(measurements);
	% now search for tag
	hasTag = false(nMeasurements,1);

	for ii = 1:nMeasurements
		mtags = ovation.asarray(measurements(ii).getAllTags());
		nTags = numel(mtags);
		if nTags > 0
			isTag = false(nTags,1);
			for jj = 1:nTags
				isTag(jj) = strcmp(tag, mtags(jj));
			end
			hasTag(ii) = any(isTag);
		end
	end
	if any(hasTag)
		measurement = measurements(hasTag);
	else
		measurement = [];
	end

	return

end

