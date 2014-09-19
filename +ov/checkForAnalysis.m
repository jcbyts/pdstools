function analysisRecord = checkForAnalysis(entity, Name, Param)
% check for analysis by protocol and protocol parameters instead of by name
% UNDER DEVELOPMENT
anRecs = ov.getAnalysis(entity);

anRecs(1).getName

Name = 'Population PSTH';
arName = [Name '_' date];

pr = context.getOrInsertProtocol(Name, '');
protocol = pr.get;

nAnRecs = numel(anRecs);

for 
singleNeuronPhys = ov.getAnalysis(entity, 'Single Neuron Physiology');

ar = entity.addAnalysisRecord(arName,...
        anRecs(1),...
        protocol,... % no protocol
        ovation.struct2map(Param));