function [bx,by, bcenters, count]=projectedHistogram(xy, n, v)
% plot the difference histogram on the unity line
% haven't figured out which options to include yet
% [bx, by, bcenters, count]=plotUnityHistogram(xy, nBins, v)
% Inputs:
%   xy [n x 2] matrix of xy points
%   v  [1 x 2] slope of line
% Outpus:
%   bx      - vector of x coordinates for plotting histogram manually
%   by      - vector of y coordinates for plotting histogram manually
%   bcenters [n x 1] bin centers
%   count    [n x 1] count at each bin
% jly wrote it

if ~exist('n', 'var') || isempty(n)
    n=20;
end

if ~exist('v','var') || isempty(v)
    v=[1 -1]; % assume it's a unity plot
end

v=v/norm(v); % line orthogonal to unity
xproj=xy*v'; % project xy points onto that line

[cnt,bins]=hist(xproj, n); % histogram the projection
count=cnt;
bcenters=bins;
cnt=cnt/sum(cnt)*mean(range(xy))/1; % normalize

off=mean(max(xy)); % offset by the xlim

% make stairs out of histogram output
[bins, cnt]=stairs(bins-mean(diff(bins))/2, cnt);
% pad appropriately
bins=[bins(1); bins; bins(end)+mean(diff(bins))*2; bins(end)+mean(diff(bins))*2; bins(1)];
cnt=[0; cnt; cnt(end); 0; 0];
% transform coordinates
bx=bins*v(1)+off+cnt;
by=bins*v(2)+off+cnt;