function experiment = getExperiment(context, projname, exname)
% GET experiment object if it exists
% example calls:
% 	1) empty = getExperiment(context)
% 		prints project and experiment names to screen, returns empty
% 	2) project = getExperiment(context, projname)
%       prints experiments for <projname> and returns project object
% 		if it exists. Otherwise it returns empty
% 	3) experiment = getExperiment(context, projname, exname)
% 		returns experiment object if it exists in projname
% 		returns project if it doesn't and project exists
% 		returns empty if neither exist
% 	4) experiment = getExperiment(project, exname)
% 		returns experiment object if it exists in projname
% 		returns project if it doesn't and project exists
% 		returns empty if neither exist
%
% Inputs:
% 		context [java obj] - an ovation context (call ovation.NewDataContext('username'))
%		projname  [string] - name of project
% 		exname 	  [string] - name of experiment
% Output:
% 		experiment [java obj] - an ovation experiment object

% 20140613 jly 	wrote it
% 20140821 jly 	package version
import ovation.*
import ov.*

if nargin < 3 || isempty(exname)
	exname = [];
	if nargin < 2 || isempty(projname)
		if isa(context, 'us.physion.ovation.domain.concrete.Project')
			projname = context;
		else
			projname = [];
		end
		if nargin < 1 || isempty(context)
			context = ovation.NewDataContext();
		end
	end
end

if isa(context, 'us.physion.ovation.domain.concrete.Project')
		startDate = exname;
		exname = projname;	
		projname = context;
end

experiment = [];


if isa(projname, 'us.physion.ovation.domain.concrete.Project')
	proj = projname;
	projname = char(proj.getName());
	nProjects = 1;
	projid 	  = 1;
else
	proj = ovation.asarray(context.getProjects());
	nProjects = numel(proj);
	projnames = cell(nProjects,1);
	for ii = 1:nProjects
		projnames{ii} = char(proj(ii).getName);
	end

	if isempty(projname)
		fprintf('Found %d Projects:\n', nProjects)
		fprintf('%s\n', projnames{:})
		projid = [];
	else
		projid = find(cellfun(@(x) strcmp(x, projname), projnames));
	end
end

if isempty(projid)
	nExperiments = zeros(nProjects,1);
	for ii = 1:nProjects
		experiments = ovation.asarray(proj(ii).getExperiments());
		nExperiments(ii) = numel(experiments);
	end
	nTotalExperiments = sum(nExperiments);
	exnames = cell(nTotalExperiments,1);
	ctr = 1;
	for ii = 1:nProjects
		experiments = ovation.asarray(proj(ii).getExperiments());
		fprintf('Project: %s has %d experiments:\n', projnames{ii}, nExperiments(ii))
		for jj = 1:nExperiments(ii)
			exnames{ctr} = char(experiments(jj).getPurpose);
			fprintf('%d)\t%s\n', jj,exnames{ctr});
			ctr = ctr + 1;
		end
	end
	
	return
else
	fprintf('found project %s.\n', projname)
	experiments = ovation.asarray(proj(projid).getExperiments());
	nExperiments = numel(experiments);
	exnames = cell(nExperiments,1);
	for ii = 1:nExperiments
		exnames{ii} = char(experiments(ii).getPurpose);
	end
end


if isempty(exname) && ~isempty(projname)
	fprintf('Project: %s has %d experiments:\n', projname, nExperiments)
	for ii = 1:nExperiments
		fprintf('%d)\t%s\n', ii,exnames{ii})
	end
	experiment = proj(projid);
	return
else
	exid = find(cellfun(@(x) strcmp(x, exname), exnames));
	if isempty(exid)
		experiment = proj(projid);
		fprintf('no experiments found matching name: %s\n', exname)
		if exist('startDate', 'var') && isa(startDate, 'org.joda.time.DateTime')
			fprintf(['\n\n*****************************************************************\n' ...
			'*****************************************************************\n\n\n']);
			fprintf('startDate object was passed in...\nAdding experiment %s now\n', exname)
			project = proj(projid);
			experiment = project.insertExperiment(exname, startDate);
		end
		return
	end
	% add experimental protocol parameters
	experiment = experiments(exid);
end
