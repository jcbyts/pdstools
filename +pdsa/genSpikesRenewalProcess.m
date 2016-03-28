function spks = genSpikesRenewalProcess(lambda, k,binsize)
% generate spike counts from a gamma renewal process
% spks = genSpikesRenewalProcess(lambda, k, dt)
% Inputs:
%   lambda  [n x 1] rate (in spikes per second)
%   k       [1 x 1] shape parameter (k=1 poisson,k < 1 more variable, k>1 more regular)
%   binsize [1 x 1] bin size (in seconds) normalizes lambda
% Output
%   spks    [n x 1] spike output
% 2015 jly wrote it

if ~exist('k', 'var')
    k = 1; % for poisson process
end

if ~exist('binsize', 'var')
    binsize = 1e-3; % assume ms bins
end
s=size(lambda);
lambda=lambda(:); % reshape into column vector
% % keep track of zeros and remove them because they mess up the mono
% l0=zeros(s);
% zeroindex=l0==0;
% lambda(zeroindex)=[]; remove
nbins       = numel(lambda);
duration    = nbins*binsize; % in seconds
time        = linspace(0,duration, nbins);

% Time rescaling gamma
intensity = cumsum(lambda)*binsize;

% monotonicity required for interp
% nonmonotonic=[0; diff(intensity)==0];

isi       = gamrnd(k,1/k,[nbins 1]);

sisi      = cumsum(isi);
sisi      = sisi(sisi<max(intensity));
y         = interp1(intensity, time, sisi);
% y         = interp1(intensity(~nonmonotonic), time(~nonmonotonic), sisi);
spks      = histc(y, time);
spks=reshape(spks, s);
