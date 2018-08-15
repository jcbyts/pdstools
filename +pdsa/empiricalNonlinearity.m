function [sprate, binCenters, praw, pspk]=empiricalNonlinearity(srate,xproj,nBins)
% check spike nonlinearity
% [sprate, binCenters, praw, pspk]=empiricalNonlinearity(srate,xproj,nBins)
%  srate [n x 1] - binned spike rate (binned spikes / binSize)
%  xproj [n x 1] - generator signal (X*w)
%  nBins [1 x 1] - resolution to plot
%
% Example:
%
% n = 1e3;
% X = randn(n,1);
% Y = poissrnd(exp(X));
% [sprate, binCenters, praw, pspk] = pdsa.empiricalNonlinearity(Y, X, 100);
% subplot(1,2,1)
% plot(binCenters, sprate, 'k');
% subplot(1,2,2)
% plot(binCenters, pspk) % only accurate if bin size <=1 ms (i.e., only 1
%                          or 0 spike possible)


 
if ~exist('nBins', 'var')
    nBins    = 50;
end

if numel(nBins) == 1 % nBins is the number of bins 
    binEdges   = (quantile(xproj, linspace(0,1,nBins+1)));
    % binCenters = (quantile(xproj, linspace(1/(2*nBins),1-(1/(2*nBins)), nBins)));
    binEdges(diff(binEdges)<=0)=[];
    binCenters=binEdges(1:end-1)+diff(binEdges)/2;
else % nBins is a vector of binCenters
    binCenters = nBins(:)'; % row vector of binCenters
    binSize = mean(diff(binCenters));
    binEdges = [binCenters - binSize/2 binCenters(end) + binSize/2];
end

% histogram xproj and kTxspk
[praw, id] = histc(xproj, binEdges);
praw = praw(1:numel(binCenters)); %praw(end) = [];

iix = praw~=0;
% binCenters(praw==0) = [];
 
nBins = numel(binCenters);
ids = 1:nBins;
sprate = zeros(nBins,1);
% ids = unique(id);
% sprate = zeros(numel(ids)-1,1);
 
for k = 1:nBins
    if praw(k) == 0
        sprate(k) = nan;
    else
        sprate(k) = mean(srate(id==ids(k)));
    end
end
 
nb=numel(binCenters);
ns=numel(sprate);
np=numel(praw);
if (nb ~= ns) || (nb ~= np)
    mn=min([nb ns np]);
    binCenters=binCenters(1:mn);
    sprate=sprate(1:mn);
    praw=praw(1:mn);
end
     
pspk = sprate./praw;