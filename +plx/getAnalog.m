function [an_info, an_data] = getAnalog(pl, analogChannels)
% Get analog data from .plx file for specific channels
% [an_info, an_data] = getAnalog(plstruct, analogChannels)
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

if ~isfield(pl.ContinuousChannels(analogChannels(1)), 'Values')
    disp('no analog data in this file')
    return
end

an_info.adfreq      = pl.ContinuousChannels(analogChannels(1)).ADFrequency;

nChannels = numel(analogChannels);
fragSize = zeros(nChannels,1);
% sum all fragments for each channels to get the total number of
% samples per channel
for ii = 1:nChannels
    fragSize(ii) = sum(pl.ContinuousChannels(analogChannels(ii)).Fragments);
end
% check for empty channels
bad = fragSize == 0;
% remove them
analogChannels(bad) = [];
minSamples = min(fragSize(~bad));
nChannels = numel(analogChannels);

an_info.nsamples    = sum(pl.ContinuousChannels(analogChannels(1)).Fragments);
an_info.timestamps  = double(pl.ContinuousChannels(analogChannels(1)).Timestamps)/pl.ADFrequency; %pl.ContinuousChannels(analogChannels(1)).ADFrequency;
an_info.fragsamples = double(pl.ContinuousChannels(analogChannels(1)).Fragments);
an_info.channels    = [pl.ContinuousChannels(analogChannels).Channel]+1; % convert from zero base
an_info.adgain      = [pl.ContinuousChannels(analogChannels).ADGain];
an_info.preampgain  = [pl.ContinuousChannels(analogChannels).PreAmpGain];


fprintf('Converting analog data to millivolts\n')

% initialize data to double (we need it in this format later)
an_data = zeros(minSamples, nChannels, 'double');
for ii = 1:nChannels
    fprintf('Channel %d of %d\n', ii, nChannels)
    thisChannel = analogChannels(ii);
    % convert values to milliVolts
    an_data(:,ii) = (pl.ContMaxMagnitudeMV*double(pl.ContinuousChannels(thisChannel).Values(1:minSamples))) ...
        / (2^(pl.BitsPerContSample-1) * pl.ContinuousChannels(thisChannel).PreAmpGain * pl.ContinuousChannels(thisChannel).ADGain);
end

fprintf('Done\n')

