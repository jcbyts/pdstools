function [spikeStructs, meanWaveforms] = concatenateSpikes(varargin)
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


%---------------------------------------------------------------------%
% plot mean waveforms
figure(100); clf
set(gca, 'Color', 'w')
for ii = 1:nSpikeStructs
	fprintf('spike %d:\t %d neurons\r')
	% figure(100+ii); clf
	subplot(1,nSpikeStructs, ii)
	plot(bsxfun(@plus, 1:n(ii), meanWaveforms{ii}')); hold on
	plot([(1:max(n-1))*nSamples(1); (1:max(n-1))*nSamples(1)], [0 max(n)+1], 'k:');

	ylim([0 max(n)+1])
	title(sprintf('spikes %d', ii))
	set(gca, 'Xtick', (	((1:max(n-1)) - 1)*nSamples(1) + 1) + nSamples(1)/2, ...
		'XtickLabel', channels-first_continuous_channel(1))
end
