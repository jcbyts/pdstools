function analysisRecord = checkForAnalysis(entity, protocolName, Param)
% check for analysis by protocol and protocol parameters instead of by name
% UNDER DEVELOPMENT
% anRecs = ov.getAnalysis(entity);
% anRecs(1).getName

import ovation.*

analysisRecord = [];

context = entity.getDataContext();
entity = context.refreshEntity(entity);

projs = ovation.asarray(entity.getProjects);
project = projs(1);

assert(numel(project) == 1, 'only one project works here')

pr = context.getOrInsertProtocol(protocolName, [protocolName '.m']);

if pr.isNew
    fprintf('New Protocol Created\r')
end

% build query options
paramFields = fieldnames(Param);
nParams = numel(paramFields);
values = cell(2,nParams);
for ii = 1:nParams
    values{1,ii} = paramFields{ii};
    tmp = Param.(paramFields{ii});
    switch class(tmp)
        case 'char'
            values{2,ii} = tmp;
        case 'double'
            values{2,ii} = num2str(tmp);
    end
end

searchStr = sprintf('param:%s=%s && ', values{:});
searchStr(end-3:end) = [];

% Search for all AnalysisRecords with parameters (note boolean logic)
results = ovation.query(context, searchStr);
nResults = numel(results);

% Iterate the results and check each record for protocol name and project
for ii = 1:nResults
    obj = results(ii);
    if (~isempty(obj) && strcmp(char(obj.getClass().getSimpleName()), 'AnalysisRecord'))
        prot = obj.getProtocol();
        if(~isempty(prot) && strcmp(char(prot.getName()), protocolName))
            if(obj.getParent().equals(project))
                disp('Found it!');
                analysisRecord = obj;
                return
            end
        end
    end
end

if isempty(analysisRecord)
    makeNew = input('do you want to add the analysis record?');

    if makeNew
        arName   = [protocolName '_' date];
        protocol = pr.get;
        an = ov.getAnalysis(project,1);
        output = an.getOutputs;
        
        
        analysisRecord = project.addAnalysisRecord(arName,...
        output,...
        protocol,...
        ovation.struct2map(Param));
    else
        return
    end
else
    return
end