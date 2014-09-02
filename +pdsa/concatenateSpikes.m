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

spikeStructs = varargin(:); % leave while debugging
nSpikeStructs = numel(spikeStructs);

n = zeros(nSpikeStructs,1);
nSamples = zeros(nSpikeStructs,1);
list = cell(nSpikeStructs,1);
channels = []; 
first_continuous_channel = zeros(nSpikeStructs,1);
remove = false(nSpikeStructs,1);
for ii = 1:nSpikeStructs
    if ~isempty(spikeStructs{ii})
	list{ii} = unique(spikeStructs{ii}.id);
	n(ii) = numel(list{ii});
	nSamples(ii) = size(spikeStructs{ii}.waveform,2);
	channels = [channels; spikeStructs{ii}.channel(:)];
    if ~isfield(spikeStructs{ii}, 'first_continuous_channel')
        spikeStructs{ii}.first_continuous_channel = 64;
    end
	first_continuous_channel(ii) = spikeStructs{ii}.first_continuous_channel;
    else
        remove(ii) = true;
    end
end

spikeStructs(remove) = [];
list(remove) = [];
n(remove) = [];
nSamples(remove) = [];
first_continuous_channel(remove) = [];
nSpikeStructs = nSpikeStructs-sum(remove);

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
		fprintf('\tneuron %d, id %d, ch: %d, snr %2.2f\r', jj, list{ii}(jj), spikeStructs{ii}.channel(jj), spikeStructs{ii}.snr(jj))
	end
end

% normalize waveform amplitude and concatenate
assert(all(nSamples==nSamples(1)), 'all waveforms must be the same size')


%% check all 3 way combos
[s1, s2, s3] = ndgrid(list{1:3});
combos = [s1(:) s2(:) s3(:)];
nCombos = size(combos,1);

M = cell(3,1);
M{1} = zeros(n(1), n(2));
M{2} = zeros(n(1), n(3));
M{3} = zeros(n(2), n(3));
% miniB = 10;
D = zeros(nCombos, 3); % stupid distance metric
for ii = 1:nCombos
    C = corrcoef([meanWaveforms{1}(combos(ii,1),:)' meanWaveforms{2}(combos(ii,2),:)' meanWaveforms{3}(combos(ii,3),:)']);
    D(ii,:) = C([2 3 6]);
end

idx = sum(D,2) < 1e-3;
combos(idx,:) = [];
D(idx,:) = [];

ind = sub2ind(size(M{1}), combos(:,1), combos(:,2));
M{1}(ind) = D(:,1);
[~, id1] = max(M{1},[], 1);
[~, id2] = max(M{1},[], 2);
nb = numel((id1));
mm = false(nb,1);
for bb = 1:nb
    mm(bb) = bb==(id2(id1(bb)));
end

nb = numel((id2));
nn = false(nb,1);
for bb = 1:nb
    nn(bb) = bb==(id1(id2(bb)));
end

[find(mm) id1(mm)']
[find(nn) id2(nn)]
%%


%---------------------------------------------------------------------%
% plot mean waveforms
figure(100); clf
set(gca, 'Color', 'w')
for ii = 1:nSpikeStructs
	fprintf('spike %d:\t %d neurons\r', ii, n(ii))
	% figure(100+ii); clf
	subplot(1,nSpikeStructs, ii)
	% mean waveforms shifted by unit #
	plot(bsxfun(@plus, 1:n(ii), meanWaveforms{ii}')); hold on
	plot([(1:max(nChannels-1))*nSamples(1); (1:max(nChannels-1))*nSamples(1)], [0 max(n)+1], 'k:');

	ylim([0 max(n)+1])
	title(sprintf('spikes %d', ii))
	set(gca, 'Xtick', (	((1:max(n-1)) - 1)*nSamples(1) + 1) + nSamples(1)/2, ...
		'XtickLabel', channels-first_continuous_channel(1))
end

fprintf('\r\r\r\r')
% for pairs at a time... will functionize 
%---------------------------------------------------------------------%
% % Waveform matching: from Tolias et al 2007
% baseStruct = 1; % start with first struct
% remaining = setdiff(1:nSpikeStructs, baseStruct);
% 
% pairwiseCombos = nchoosek(1:nSpikeStructs,2);
% % only keep linking pairs
% [~,keepIndex] = unique(pairwiseCombos(:,1));
% pairwiseCombos = pairwiseCombos(keepIndex,:);
% nPairs = size(pairwiseCombos,1);
% for ii = 1:nPairs
% 	[a, d1, d2] = pdsa.waveformMatch(meanWaveforms{pairwiseCombos(ii,1)}', meanWaveforms{pairwiseCombos(ii,2)}');
% 	[~, id] = sort(d2-d1, 2);
% 	[un, unidx] = unique(id(:,end));
% 	if numel(unique(id(:,end))) == numel(id(:,end))
% 		unitMatch{ii} = id(:,end);
% 	else
% 		commandwindow
% 		fprintf('there are duplicate units')
% 		unitMatch{ii} = nan(n(pairwiseCombos(ii,1)),1);
% 		for jj = 1:n(pairwiseCombos(ii,1))
% 			tmp = input(sprintf('unit %d spikes %d is which in spikes %d', jj, pairwiseCombos(ii,1), pairwiseCombos(ii,2)));
% 			if ~isempty(tmp)
% 				unitMatch{ii}(jj) = tmp;
% 			end
% 		end
% 	end
% end

% keyboard

unitMatch = cell(nSpikeStructs,1);
for ii = 1:nSpikeStructs
    for jj = 1:n(ii)
    	figure(100)
        tmp = input(sprintf('unit %d spikes %d newlabel', jj, ii));
        if ~isempty(tmp)
            unitMatch{ii}(jj) = tmp;
        end
    end
end
for ii = 1:nSpikeStructs
	tmpSpikeStruct = spikeStructs{ii};
	spikeStructs{ii}.snr 	 = nan(max(unitMatch{ii}),1);
	spikeStructs{ii}.channel = nan(max(unitMatch{ii}),1);
	for jj = 1:n(ii)
        if unitMatch{ii}(jj)~=0
            spikeStructs{ii}.id(tmpSpikeStruct.id==jj) = unitMatch{ii}(jj);
            spikeStructs{ii}.snr(unitMatch{ii}(jj)) = tmpSpikeStruct.snr(jj);
            spikeStructs{ii}.channel(unitMatch{ii}(jj)) = tmpSpikeStruct.channel(jj);
        end
	end
end

spikesNew = [];
% % idea is simple: get back an a, d1 and d2;  Loop over units in spikes1 and check if any of them exist in spikes 2.
% % If they do, great, combine them. If they don't, add them as new units.
% % Step 1: get start and stop time of each spikes struct -- the first one is zero to last spike.
% tStart = zeros(nSpikeStructs,1);
% tStop  = zeros(nSpikeStructs,1);

% tStop(1) = spikeStructs{1}.time(end)+10; % add 10 seconds to first file
% for ii = 2:nSpikeStructs
% 	tStart(ii) = tStop(ii-1);
% 	tStop(ii)  = tStart(ii) + spikeStructs{ii}.time(end) + 10;
% end

% % Loop over units in spikes1 and find best match in spikes 2
% unitMatch = zeros(n(1),3); %[id a(id) d1(id) d2(id)]
% for ii = 1:n(1)
% 		[~, id1] = min(d1(ii,:));
% 		[~, id2] = max(d2(ii,:));
% 		if id1 == id2
% 			unitMatch(ii,:) = [id1 d1(ii,id1) d2(ii,id2)];
% 		end

% end

% % Check matches and correct for any duplicates
% uniqueMatches = unique(unitMatch(:,1));
% for kk = 1:numel(uniqueMatches)
% 	unitCheck = unitMatch(:,1) == uniqueMatches(kk);
% 	if any(unitCheck)
% 		[~,idx1] = min(unitMatch(:,2));
% 		[~,idx2] = max(unitMatch(:,3));	
% 		if idx1 == idx1
% 			xid = setdiff(find(unitCheck), idx1);
% 			unitMatch(xid,:) = 0;
% 		end
% 	end
% end

% spikesTmp1 = spikeStructs{1};
% spikesTmp2 = spikeStructs{2};
% spikesTmp2.id = spikeStructs{2}.id + n(1);
% spikesTmp2.time = spikesTmp2.time + tStart(2);
% for ii = 1:n(1)
% 	if unitMatch(ii,1) ~=0
% 		idx = spikesTmp2.id == unitMatch(ii,1) + n(1);
% 		spikesTmp2.id(idx) = ii;
% 		spikesTmp2.channel(unitMatch(ii,1)) = 0;
% 		spikesTmp2.snr(unitMatch(ii,1)) = 0;
% 	end
% end



% spikesNew = struct('time', [spikesTmp1.time; spikesTmp2.time], ...
% 	'id', [spikesTmp1.id; spikesTmp2.id], ...
% 	'waveform', [spikesTmp1.waveform; spikesTmp2.waveform], ...
% 	'channel', [spikesTmp1.channel spikesTmp2.channel], ...
% 	'snr', [spikesTmp1.snr spikesTmp2.snr]);

% spikesNew.channel(spikesNew.channel == 0) = [];
% spikesNew.snr(spikesNew.snr == 0) = [];

% spikesTmp2.time = spikesTmp2.time - tStart(2);
% spikeStructs{1} = spikesTmp1;
% spikeStructs{2} = spikesTmp2;

% % spikeStructs{end+1} = spikesNew;
% % nSpikeStructs = numel(spikeStructs);



% %---------------------------------------------------------------------%
% % plot mean waveforms
% figure(100); clf
% set(gca, 'Color', 'w')
% for ii = 1:nSpikeStructs
% 	fprintf('spike %d:\t %d neurons\r')
% 	% figure(100+ii); clf
% 	subplot(1,nSpikeStructs, ii)
% 	% mean waveforms shifted by unit #
% 	plot(bsxfun(@plus, 1:n(ii), meanWaveforms{ii}')); hold on
% 	plot([(1:max(n-1))*nSamples(1); (1:max(n-1))*nSamples(1)], [0 max(n)+1], 'k:');

% 	ylim([0 max(n)+1])
% 	title(sprintf('spikes %d', ii))
% 	set(gca, 'Xtick', (	((1:max(n-1)) - 1)*nSamples(1) + 1) + nSamples(1)/2, ...
% 		'XtickLabel', channels-first_continuous_channel(1))
% end
