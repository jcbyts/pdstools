function [spikeStructs, spikeStructsNew] = concatenateSpikes(varargin)
% Concatenate two spike time structures
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
spikeStructsNew = spikeStructs;
nSpikeStructs = numel(spikeStructs);
wf = zeros(nSpikeStructs,1);
for ii = 1:nSpikeStructs
    wf(ii) = size(spikeStructs{ii}.waveform,2);
end

if ~all(wf==wf(1))
    disp('mistmatched waveform samples')
    uwf = unique(wf);
    assert(numel(uwf)==2, 'wtf. only two sampling rates')
    nwf = lcm(uwf(1),uwf(2));
    
    for ii = 1:nSpikeStructs
        spikeStructs{ii}.waveform_old = spikeStructs{ii}.waveform;
        nw = size(spikeStructs{ii}.waveform,1);
        nws = size(spikeStructs{ii}.waveform,2);
        spikeStructs{ii}.waveform = zeros(nw, nwf);
        for jj = 1:nw
            spikeStructs{ii}.waveform(jj,:) = interp(spikeStructs{ii}.waveform_old(jj,:), nwf/nws);
        end
        fprintf('%2.0f done\n', ii/nSpikeStructs*100)
    end
end

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
        channels = [channels; spikeStructs{ii}.channel(list{ii})']; %#ok<*AGROW>
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
meanWaveforms = cell(nSpikeStructs,1);
for ii = 1:nSpikeStructs
	meanWaveforms{ii} = zeros(n(ii), nChannels*nSamples(ii));
	for jj = 1:n(ii)
		ch = spikeStructs{ii}.channel(list{ii}(jj)) - spikeStructs{ii}.first_continuous_channel;
        if ~isnan(ch)
            ch = find(channels-first_continuous_channel(ii)==ch);
            idx = ((ch-1)*nSamples(ii) + 1):(ch*nSamples(ii));
            mw=mean(spikeStructs{ii}.waveform(spikeStructs{ii}.id==list{ii}(jj),:));
            mw = mw - mean(mw);
            mw = mw/norm(mw);
            meanWaveforms{ii}(jj,idx) = mw;
        end
	end
end


fprintf('************************************************\n')
fprintf('found %d spike files with %d neurons\n', nSpikeStructs, max(n))
for ii = 1:nSpikeStructs
	fprintf('spike%d:\n', ii)
	for jj = 1:n(ii)
		fprintf('\tneuron %02.0f, id %02.0f, ch: %d, snr %2.2f\n', jj, list{ii}(jj), spikeStructs{ii}.channel(jj), spikeStructs{ii}.snr(jj))
	end
end


%---------------------------------------------------------------------%
% plot mean waveforms

figure(100); clf
set(gca, 'Color', 'w')
for ii = 1:nSpikeStructs
	fprintf('spike %d:\t %d neurons\n', ii, n(ii))
	% figure(100+ii); clf
	subplot(1,nSpikeStructs, ii)
    
% 	plot(bsxfun(@plus,1:n(ii), meanWaveforms{ii}(1:n(ii),:)')); hold on
    plot(bsxfun(@plus,list{ii}', meanWaveforms{ii}(1:n(ii),:)')); hold on
	plot([(1:max(nChannels))*nSamples(1); (1:max(nChannels))*nSamples(ii)], [0 max(n)+1], 'k:');
    
% 	ylim([0 max(n)+1])
    ylim([0 max(cell2mat(list(:)))+1])
    
	title(sprintf('spikes %d', ii))
%     channels(~isnan(channels))-first_continuous_channel(ii)
    xax = (	((1:(max(n)))-1)*nSamples(1)) + nSamples(1)/2;
	set(gca, 'Xtick', xax(1:end), ...
		'XtickLabel', channels-first_continuous_channel(ii))
%     xlim([1 (max(channels) - first_continuous_channel(ii))*nSamples(ii)]) 
    xlim([1 size(meanWaveforms{ii},2)])
end

continueMatch = input('try waveform match? (1 or 0)');
if ~continueMatch
    return
end
% normalize waveform amplitude and concatenate
assert(all(nSamples==nSamples(1)), 'all waveforms must be the same size')


%% check all 3 way combos
units = cellfun(@(x) 1:x, num2cell(n), 'Uniformoutput', false);
switch nSpikeStructs
    case 2
%         [s1, s2] = ndgrid(list{1:nSpikeStructs});
        [s1, s2] = ndgrid(units{:});
        combos = [s1(:) s2(:)];
    case 3
%         [s1, s2, s3] = ndgrid(list{1:nSpikeStructs});
        [s1, s2, s3] = ndgrid(units{:});
        combos = [s1(:) s2(:) s3(:)];
end
nCombos = size(combos,1);

pairwise = nchoosek(1:nSpikeStructs,2);
nPairs   = size(pairwise,1);
M = cell(nPairs,1);
for ii = 1:nPairs
	M{ii} = zeros(n(pairwise(ii,1)), n(pairwise(ii,2)));
end

D = zeros(nCombos, nPairs); % stupid distance metric
for ii = 1:nCombos
    for jj = 1:nPairs
        [~,~,d2] = pdsa.waveformMatch(meanWaveforms{pairwise(jj,1)}(combos(ii,pairwise(jj,1)),:)', meanWaveforms{pairwise(jj,2)}(combos(ii,pairwise(jj,2)),:)');
        D(ii,jj) = d2;
    end
end

% remove all near zeros
idx = sum(D,2) < 1e-3;
combos(idx,:) = [];
D(idx,:) = [];
match = [];
%%
for ii = 1:nPairs
    ind = sub2ind(size(M{ii}), combos(:,pairwise(ii,1)), combos(:,pairwise(ii,2)));
    M{ii}(ind) = D(:,ii);
    [~, id2] = max(M{ii},[], 1);
    [~, id1] = max(M{ii},[], 2);
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
    
    
    
    if isempty(match)
        match = [match; [find(mm) id1(mm) nan(numel(find(mm)),1)]];
    else
        [i,~] = find(bsxfun(@eq, match(:,pairwise(ii,1)),find(mm)'));
        match(i, pairwise(ii,2)) = id1(mm);
    end
    
    
    %%
    [i,~] = find(bsxfun(@eq, match(:,pairwise(ii,2)), find(nn)'));
    assert(all(unique(match(i,pairwise(ii,2)))==find(nn)), 'they didn''t match')
    
    % find neurons from file 1 that don't match
    unmatched1 = setdiff(1:n(pairwise(ii,1)), match(:,pairwise(ii,1)));
    unmatched2 = setdiff(1:n(pairwise(ii,2)), match(:,pairwise(ii,2)));
    tmp = nan(numel(unmatched1)+numel(unmatched2),3);
    tmp(1:numel(unmatched1),pairwise(ii,1)) = unmatched1;
    tmp(numel(unmatched1)+1:end,pairwise(ii,2)) = unmatched2;
    match = [match; tmp];
end

%%
%---------------------------------------------------------------------%
% plot mean waveforms
figure(100); clf
set(gca, 'Color', 'w')
for ii = 1:nSpikeStructs
	fprintf('spike %d:\t %d neurons\n', ii, n(ii))
	% figure(100+ii); clf
	subplot(1,nSpikeStructs, ii)
	% mean waveforms shifted by unit #
    tind = ~isnan(match(:,ii));
    
	plot(bsxfun(@plus,find(tind)', meanWaveforms{ii}(match(tind,ii),:)')); hold on
	plot([(1:max(nChannels))*nSamples(1); (1:max(nChannels))*nSamples(ii)], [0 max(size(match,1))+1], 'k:');

	ylim([0 size(match,1)+1])
    
	title(sprintf('spikes %d', ii))
    xax = (	((1:(max(n)))-1)*nSamples(1)) + nSamples(1)/2;
%     xax = (	((1:(max(n)-1)) - 1)*nSamples(1)) + nSamples(1)/2;
	set(gca, 'Xtick', xax(1:end), ...
		'XtickLabel', channels-first_continuous_channel(ii))
    xlim([1 size(meanWaveforms{ii},2)]) 
%     xlim([1 xax(end)])
end

nNeurons = size(match,1);
spikesNew = repmat(struct('time', [], 'id', [], 'waveform', [], 'channel', nan(1,nNeurons), 'snr', nan(1,nNeurons), 'first_continuous_channel', 64, 'continuous_only', 1), nSpikeStructs,1);

for ii = 1:nSpikeStructs
%     chans = spikeStructs{ii}.channel(~isnan(spikeStructs{ii}.channel));
%     snr   = spikeStructs{ii}.snr(~isnan(spikeStructs{ii}.snr));
    for jj = 1:nNeurons
        if ~isnan(match(jj,ii))
            
            ind = spikeStructs{ii}.id==list{ii}(match(jj,ii));
            spks = spikeStructs{ii}.time(ind);
            spikesNew(ii).time = [spikesNew(ii).time; spks(:)];
            spikesNew(ii).id   = [spikesNew(ii).id; ones(numel(spks),1)*jj];
            spikesNew(ii).waveform = [spikesNew(ii).waveform; spikeStructs{ii}.waveform(ind,:)];
            spikesNew(ii).channel(jj) = spikeStructs{ii}.channel(list{ii}(match(jj,ii)));
%             spikesNew(ii).channel(jj) = chans(match(jj,ii));
            spikesNew(ii).snr(jj) = spikeStructs{ii}.snr(list{ii}(match(jj,ii)));
%             spikesNew(ii).snr(jj) = snr(match(jj,ii));
        end
        if isfield(spikeStructs{ii}, 'first_continuous_channel')
            spikesNew(ii).first_continuous_channel = spikeStructs{ii}.first_continuous_channel;
        end
    end
end
%%
for ii = 1:nSpikeStructs
    spikeStructsNew{ii} = spikesNew(ii);
end
fprintf('\n\n\n\n')

% % good = input('is it good?');
% % if good
% %     
% %     
% %     
% %     
% %     
% %     
% %     
% %     
% %     
% %     
% %     
% %     
% % end    
%     
%     
%     
%     
%     
%     
% %%
% % for pairs at a time... will functionize 
% %---------------------------------------------------------------------%
% % % Waveform matching: from Tolias et al 2007
% % baseStruct = 1; % start with first struct
% % remaining = setdiff(1:nSpikeStructs, baseStruct);
% % 
% % pairwiseCombos = nchoosek(1:nSpikeStructs,2);
% % % only keep linking pairs
% % [~,keepIndex] = unique(pairwiseCombos(:,1));
% % pairwiseCombos = pairwiseCombos(keepIndex,:);
% % nPairs = size(pairwiseCombos,1);
% % for ii = 1:nPairs
% % 	[a, d1, d2] = pdsa.waveformMatch(meanWaveforms{pairwiseCombos(ii,1)}', meanWaveforms{pairwiseCombos(ii,2)}');
% % 	[~, id] = sort(d2-d1, 2);
% % 	[un, unidx] = unique(id(:,end));
% % 	if numel(unique(id(:,end))) == numel(id(:,end))
% % 		unitMatch{ii} = id(:,end);
% % 	else
% % 		commandwindow
% % 		fprintf('there are duplicate units')
% % 		unitMatch{ii} = nan(n(pairwiseCombos(ii,1)),1);
% % 		for jj = 1:n(pairwiseCombos(ii,1))
% % 			tmp = input(sprintf('unit %d spikes %d is which in spikes %d', jj, pairwiseCombos(ii,1), pairwiseCombos(ii,2)));
% % 			if ~isempty(tmp)
% % 				unitMatch{ii}(jj) = tmp;
% % 			end
% % 		end
% % 	end
% % end
% 
% % keyboard
% 
% unitMatch = cell(nSpikeStructs,1);
% for ii = 1:nSpikeStructs
%     for jj = 1:n(ii)
%     	figure(100)
%         tmp = input(sprintf('unit %d spikes %d newlabel', jj, ii));
%         if ~isempty(tmp)
%             unitMatch{ii}(jj) = tmp;
%         end
%     end
% end
% for ii = 1:nSpikeStructs
% 	tmpSpikeStruct = spikeStructs{ii};
% 	spikeStructs{ii}.snr 	 = nan(max(unitMatch{ii}),1);
% 	spikeStructs{ii}.channel = nan(max(unitMatch{ii}),1);
% 	for jj = 1:n(ii)
%         if unitMatch{ii}(jj)~=0
%             spikeStructs{ii}.id(tmpSpikeStruct.id==jj) = unitMatch{ii}(jj);
%             spikeStructs{ii}.snr(unitMatch{ii}(jj)) = tmpSpikeStruct.snr(jj);
%             spikeStructs{ii}.channel(unitMatch{ii}(jj)) = tmpSpikeStruct.channel(jj);
%         end
% 	end
% end
% 
% spikesNew = [];
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
% 	fprintf('spike %d:\t %d neurons\n')
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
