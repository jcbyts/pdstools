function [bx,by, bcenters, count, offset]=projectedHistogram(xy, n, v, offset, normzer, rnge)
% plot the difference histogram on the unity line
% haven't figured out which options to include yet
% [bx, by, bcenters, count]=plotUnityHistogram(xy, nBins, v, offset, normzer, rnge)
% Inputs:
%   xy      [n x 2] matrix of xy points
%   n       [1 x 1] number of bins (or vector of bin edges)
%   v       [1 x 2] slope of line
%   offset  [1 x 2] xy offset for histogram
%   normzer [1 x 1] normalization for count. We have to scale this thing
%   rnge    [1 x 1] scale for yaxis of histogram
% Outpus:
%   bx      - vector of x coordinates for plotting histogram manually
%   by      - vector of y coordinates for plotting histogram manually
%   bcenters [n x 1] bin centers
%   count    [n x 1] count at each bin
%   offset   [1 x 2] histogram 0,0 coordinate
% jly wrote it

if ~exist('n', 'var') || isempty(n)
    n=20;
end

if ~exist('v','var') || isempty(v)
    v=[1 -1]; % assume it's a unity plot
end

if ~exist('offset', 'var') || isempty(offset)
    offset=mean(max(xy)); % offset by the xlim
    disp(offset)
end

if ~exist('rnge', 'var') || isempty(rnge)
    rnge=mean(range(xy))/1;
end

v=v/norm(v); % line orthogonal to unity
xproj=xy*v'; % project xy points onto that line

if numel(n)>1
    n=n(:)';
    bw=diff(n)/2;
    bw=[-bw -bw(end)];
    bedges=[n+bw n(end)-bw(end)];
    bins=n;
    cnt=histc(xproj, bedges);
    cnt=cnt(1:end-1);
    cnt=cnt(:)';
else
    [cnt,bins]=hist(xproj, n); % histogram the projection
end
count=cnt;
bcenters=bins;
if ~exist('normzer', 'var') || isempty(normzer)
    normzer=sum(cnt);
end
cnt=cnt/normzer*rnge; % normalize



% make stairs out of histogram output
[bins, cnt]=stairs(bins-mean(diff(bins))/2, cnt);
% pad appropriately
bins=[bins(1); bins; bins(end)+mean(diff(bins))*2; bins(end)+mean(diff(bins))*2; bins(1)];
cnt=[0; cnt; cnt(end); 0; 0];
% transform coordinates
bx=bins*v(1)+offset+cnt;
by=bins*v(2)+offset+cnt;