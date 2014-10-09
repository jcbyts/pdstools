function param = getInfo(context, projname, exname, info)
% GET or ADD protocol parameters
% param = ovAddInfo(context, projname, exname, info)
% getInfo can ADD protocol parameters or GET protocol parameters
% Inputs: varargin
% 	ov_addInfo(experiment, [info])
% 		experiment [java obj] 
% 		info 		 [struct] (optional) - adds info struct to experiment
% 							  returns info struct if empty or absent
% 	ov_addInfo(context, projname, exname, [info])
%		context [java obj] - ovation context (ovation.NewDataContext())
% 	 	 projname [string] - project name
% 	   	   exname [string] - experiment name
% 		 	 info [struct] - fieldnames get added as protocol parameters
% Output:
% 		param [struct] - fieldnames are protocol parameters
% 
% if info is empty ov_addInfo will check for existing parameters and
% return them

% 20140613 jly 	wrote it

if nargin < 4 || ~isstruct(info)
	info = [];
	if nargin < 3 || isempty(exname)
		exname = [];
		if nargin < 2 || isempty(projname)
			projname = [];
			if nargin < 1 || isempty(context)
				context = ovation.NewDataContext();
			end
		end
	end
end

if isa(context, 'us.physion.ovation.domain.concrete.Experiment')
	experiment = context;
	if isstruct(projname)
		info = projname;
	end
	projects = ovation.asarray(experiment.getProjects());
	nProjects = numel(projects);
	projname = char(projects(1).getName());
	exname   = char(experiment.getPurpose());

elseif isa(context, 'us.physion.ovation.domain.concrete.AnalysisRecord')
	if isstruct(projname)
		info = projname;
	end
	projects = ovation.asarray(context.getProjects);
	projname = char(projects(1).getName());
	experiments = ovation.asarray(context.getExperiments);
	exname = char(experiments(1).getPurpose);
	experiment = context;
else
	experiment = ov_getExperiment(context, projname, exname);
end

if isempty(experiment)
	param = [];
	return
end

% check for existing parameters
param = ovation.map2struct(experiment.getProtocolParameters());
tags  = ovation.asarray(experiment.getAllTags());
nParam = numel(param);
nTags  = numel(tags);
infofields = fieldnames(info);
paramfields = fieldnames(param);
ismatch = false(numel(paramfields), numel(infofields));
for ii = 1:numel(paramfields)
	for jj = 1:numel(infofields)
		ismatch(ii,jj) = strcmp(paramfields{ii}, infofields{jj});
	end
end

info = rmfield(info, infofields(find(sum(ismatch))));

if nParam > 0
	params = fieldnames(param);
	if exist('info', 'var') && isstruct(info) && ~isempty(info) && (isempty(params) || ~ovation.structCmp(param, info)) 
		fields = fieldnames(info);
		nFields = numel(fields);
		for ii = 1:nFields
			if ~isempty(info.(fields{ii}))
				experiment.addProtocolParameter(fields{ii}, info.(fields{ii}));
			end
		end
		param = ovation.map2struct(experiment.getProtocolParameters());
		params = fieldnames(param);
	end
	
	nParams = numel(params);
	fprintf('Experiment: %s has %d protocol parameters\n', exname, nParams)
	for ii = 1:nParams
		fprintf('%d)%s\n',ii,params{ii})
	end
	
end

param = ovation.map2struct(experiment.getProtocolParameters());
params = fieldnames(param);
nParams = numel(params);
fprintf('Found %d protocol parameters\n', nParams)
for ii = 1:nParams
	if isa(param.(params{ii}), 'java.lang.String[]')
		nP = numel(param.(params{ii}));
		fieldData = cell(nP,1);
		for jj = 1:nP
			fieldData{jj} = char(param.(params{ii})(jj));
			tagCompare = false(nTags,1);
			for kk = 1:nTags
				tagCompare(kk) = strcmp(char(param.(params{ii})(jj)), char(tags(kk)));
			end
			if ~any(tagCompare)
				experiment.addTag(char(param.(params{ii})(jj)));
			end
		end
		param = rmfield(param, (params{ii}));
		param.(params{ii}) = fieldData;
	end
end


if isempty(params)
	fprintf('params are still empty\n')
	return
end
