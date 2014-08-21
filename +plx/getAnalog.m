function [an_info, an_data] = plx_getAnalog(plxname, channels)
%[an_info, an_data] = plx_getAnalog(plxname, channels)
% inputs:
%   plxname     [string] - full/path/to/plxname.plx
%   channels     [n x 1] - vector of channels to pull. They must all have
%                          the same sampling rate.
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

an_data = [];
an_info = [];
if nargin < 2
    help plx_getAnalog
    return
end

% get channel map to identify possible analog channels
[~, adc]    = plx_ad_chanmap(plxname);
[~, sc]     = plx_adchan_freqs(plxname); %= plx_adchan_samplecounts(plxname);

% get channels as indices
nchannels = numel(channels);
ichannels = zeros(nchannels, 1);
for ch = 1:nchannels
    ichannels(ch) = find(channels(ch)-1==adc);
end

assert(all(sc(ichannels)==sc(ichannels(1))), 'all channels specified must have the same sampling rate')

% get sampling information for base channel. All other channels should have
% the same info.
[adfreq, n, ts, fn] = plx_ad_gap_info(plxname, adc(ichannels(1)));

% get analog data in mV
startCount = 1;
endCount  = n;

% cleanup
an_info.adfreq = adfreq;
an_info.nsamples = n;
an_info.timestamps = ts;
an_info.fragsamples = fn;
an_info.channels  = channels;

if nargout > 1
    an_data = zeros(n,nchannels, 'single');
    fprintf('pulling analog data from %d channels\r', nchannels)
    for ch = 1:nchannels
        fprintf('channel %d\r', channels(ch))
        [~, ~, an_data(:,ch)] = plx_ad_span_v(plxname, channels(ch)-1, startCount, endCount);
    end
end

