function [an_info, an_data] = getAnalog(pl, analogChannels)
%[an_info, an_data] = getAnalog(plstruct, analogChannels)
% inputs:
%   plstruct 	[struct] - fullread output of readPLXc
%   channels    [double] - vector of channels to read.
%
% outputs:
%   an_data     [m x n] - (samples x channels) the data (single precision)
%   an_info    [struct] - info about sampling (so time can be
%                         reconstructed)
%       .adfreq [1 x 1] - sampling rate
%     .nsamples [1 x 1] - number of total samples
%   .timestamps [k x 1] - samples come in blocks. if the plexon recording
%                         was paused at any point, this indicates the
%                         timestamp of the first sample of each fragment.
%  .fragsamples [k x 1] - number of samples per recording fragment
%     .channels [n x 1] - vector of channel id for each column of an_data

% 20140531  jly     wrote it
% 20140820  jly     changed to readPLX version

import plx.*

an_info  = struct();
an_data  = [];
if nargin < 1
    help getAnalog
    return
end

%--------------------------------------------------------------------------------------------%
%% make analog struct: this part is tricky to make a memory efficient version



analogSamplingRates  = [pl.ContinuousChannels(analogChannels).ADFrequency]';
analogRates = unique(analogSamplingRates); % to sort into lfp and analog
% if only one rate exists, only make the lfp struct
nRates = numel(analogRates);
if nRates > 1
	fprintf('must only have one sampling rate across all the channels passed in')
	return
end

an_info.adfreq    = pl.ContinuousChannels(analogChannels(1)).ADFrequency;
an_info.nsamples = sum(pl.ContinuousChannels(analogChannels(1)).Fragments);
an_info.timestamps  = double(pl.ContinuousChannels(analogChannels(1)).Timestamps)/pl.ADFrequency; %pl.ContinuousChannels(analogChannels(1)).ADFrequency;
an_info.fragsamples = double(pl.ContinuousChannels(analogChannels(1)).Fragments);
an_info.channels    = [pl.ContinuousChannels(analogChannels).Channel];
an_info.adgain      = [pl.ContinuousChannels(analogChannels(analogChannels)).ADGain];
an_info.preampgain  = [pl.ContinuousChannels(analogChannels(analogChannels)).PreAmpGain];


try 
	an_data = [pl.ContinuousChannels(analogChannels).Values];
catch me
	fprintf('failed to concatenate. proceeding the slow way\n')
	nChannels = numel(analogChannels);
	fragSize = zeros(nChannels,1);
	for ii = 1:nChannels
		fragSize(ii) = sum(pl.ContinuousChannels(ii).Fragments);
	end
	an_data = zeros(min(fragSize), nChannels, 'int16');
	for ii = 1:nChannels
		% an_data = [pl.ContinuousChannels(analogChannels).Values];
		an_data(:,ii) = pl.ContinuousChannels(ii).Values(1:min(fragSize));
	end
end
