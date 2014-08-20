function [spikeStructs, spikesNew] = concatenateSpikes(varargin)
% concatenate two spike time structures
% example call:
% 	spikesNew = concatenateSpikes(spikes1, spikes2)
% INPUTS:
% 	spikes1 [struct]
% 		.times
% 		.id
% 		.waveform
% 		.channel
% 		.snr
% 		.first_continuous_channel

spikeStructs = varargin(1:2); % leave while debugging
nSpikeStructs = numel(spikeStructs);

n = zeros(nSpikeStructs,1);
nSamples = zeros(nSpikeStructs,1);
list = cell(nSpikeStructs,1);
channels = []; 
first_continuous_channel = zeros(nSpikeStructs,1);
for ii = 1:nSpikeStructs
	list{ii} = unique(spikeStructs{ii}.id);
	n(ii) = numel(list{ii});
	nSamples(ii) = size(spikeStructs{ii}.waveform,2);
	channels = [channels; spikeStructs{ii}.channel(:)];
	first_continuous_channel(ii) = spikeStructs{ii}.first_continuous_channel;
end


channels  = unique(channels);
nChannels = numel(channels);

for ii = 1:nSpikeStructs
	meanWaveforms{ii} = zeros(n(ii), nChannels*nSamples(1));
	for jj = 1:n(ii)
		ch = spikeStructs{ii}.channel(jj) - spikeStructs{ii}.first_continuous_channel;
		ch = find(channels-first_continuous_channel(ii)==ch);
		% idx = ((jj-1).*32 + 1):(jj.*nSamples(1));
		idx = ((ch-1).*32 + 1):(ch.*nSamples(1));
		mw=mean(spikeStructs{ii}.waveform(spikeStructs{ii}.id==list{ii}(jj),:));
		mw = mw - mean(mw);
		mw = mw/norm(mw);
		meanWaveforms{ii}(jj,idx) = mw;
	end
end


fprintf('************************************************\r')
fprintf('found %d spike files with %d neurons\r', nSpikeStructs, max(n))
for ii = 1:nSpikeStructs
	fprintf('spike%d:\r', ii)
	for jj = 1:n(ii)
		fprintf('\tneuron %d, id %d, ch: %d, snr %2.2f\r', jj, spikeStructs{ii}.id(jj), spikeStructs{ii}.channel(jj), spikeStructs{ii}.snr(jj))
	end
end

% normalize waveform amplitude and concatenate
assert(all(nSamples==nSamples(1)), 'all waveforms must be the same size')

% for pairs at a time... will functionize 
%---------------------------------------------------------------------%
% Waveform matching: from Tolias et al 2007
baseStruct = 1; % start with first struct
remaining = setdiff(1:nSpikeStructs, baseStruct);

[a, d1, d2] = pdsa.waveformMatch(meanWaveforms{1}', meanWaveforms{2}');

% idea is simple: get back an a, d1 and d2;  Loop over units in spikes1 and check if any of them exist in spikes 2.
% If they do, great, combine them. If they don't, add them as new units.
% Step 1: get start and stop time of each spikes struct -- the first one is zero to last spike.
tStart = zeros(nSpikeStructs,1);
tStop  = zeros(nSpikeStructs,1);

tStop(1) = spikeStructs{1}.time(end)+10; % add 10 seconds to first file
for ii = 2:nSpikeStructs
	tStart(ii) = tStop(ii-1);
	tStop(ii)  = tStart(ii) + spikeStructs{ii}.time(end) + 10;
end

% Loop over units in spikes1 and find best match in spikes 2
unitMatch = zeros(n(1),3); %[id a(id) d1(id) d2(id)]
for ii = 1:n(1)
		[~, id1] = min(d1(ii,:));
		[~, id2] = max(d2(ii,:));
		if id1 == id2
			unitMatch(ii,:) = [id1 d1(ii,id1) d2(ii,id2)];
		end

end

% Check matches and correct for any duplicates
uniqueMatches = unique(unitMatch(:,1));
for kk = 1:numel(uniqueMatches)
	unitCheck = unitMatch(:,1) == uniqueMatches(kk);
	if any(unitCheck)
		[~,idx1] = min(unitMatch(:,2));
		[~,idx2] = max(unitMatch(:,3));	
		if idx1 == idx1
			xid = setdiff(find(unitCheck), idx1);
			unitMatch(xid,:) = 0;
		end
	end
end

spikesTmp1 = spikeStructs{1};
spikesTmp2 = spikeStructs{2};
spikesTmp2.id = spikeStructs{2}.id + n(1);
spikesTmp2.time = spikesTmp2.time + tStart(2);
for ii = 1:n(1)
	if unitMatch(ii,1) ~=0
		idx = spikesTmp2.id == unitMatch(ii,1) + n(1);
		spikesTmp2.id(idx) = ii;
		spikesTmp2.channel(unitMatch(ii,1)) = 0;
		spikesTmp2.snr(unitMatch(ii,1)) = 0;
	end
end

spikesNew = struct('time', [spikesTmp1.time; spikesTmp2.time], ...
	'id', [spikesTmp1.id; spikesTmp2.id], ...
	'waveform', [spikesTmp1.waveform; spikesTmp2.waveform], ...
	'channel', [spikesTmp1.channel spikesTmp2.channel], ...
	'snr', [spikesTmp1.snr spikesTmp2.snr]);

spikesNew.channel(spikesNew.channel == 0) = [];
spikesNew.snr(spikesNew.snr == 0) = [];


% spikeStructs{end+1} = spikesNew;
% nSpikeStructs = numel(spikeStructs);



%---------------------------------------------------------------------%
% plot mean waveforms
figure(100); clf
set(gca, 'Color', 'w')
for ii = 1:nSpikeStructs
	fprintf('spike %d:\t %d neurons\r')
	% figure(100+ii); clf
	subplot(1,nSpikeStructs, ii)
	% mean waveforms shifted by unit #
	plot(bsxfun(@plus, 1:n(ii), meanWaveforms{ii}')); hold on
	plot([(1:max(n-1))*nSamples(1); (1:max(n-1))*nSamples(1)], [0 max(n)+1], 'k:');

	ylim([0 max(n)+1])
	title(sprintf('spikes %d', ii))
	set(gca, 'Xtick', (	((1:max(n-1)) - 1)*nSamples(1) + 1) + nSamples(1)/2, ...
		'XtickLabel', channels-first_continuous_channel(1))
end
