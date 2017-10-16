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
 
binEdges   = (quantile(xproj, linspace(0,1,nBins+1)));
% binCenters = (quantile(xproj, linspace(1/(2*nBins),1-(1/(2*nBins)), nBins)));
 
binEdges(diff(binEdges)<=0)=[];
binCenters=binEdges(1:end-1)+diff(binEdges)/2;
% binCenters(diff(binCenters)<=0)=[];
% histogram xproj and kTxspk
[praw, id] = histc(xproj, binEdges); 
praw = praw(1:numel(binCenters)); praw(end) = [];
binCenters(praw==0) = [];
 
ids = unique(id);
sprate = zeros(numel(ids)-1,1);
 
for k = 1:numel(ids)-1
    sprate(k) = mean(srate(id==ids(k)));
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