%% test scatter histogram
import pdsa.*
n1=100;
n2=50;
% make two distributions
xy1=mvnrnd([1 0], diag([1 1]), n1);
xy2=mvnrnd([.5 .5], [1 .1; .1 2], n2);

% test for two distributions
nBins=25;

xy=[xy1; xy2]; % concatenate distributions to find coordinates
[bx,~, bcenters, count, offset]=projectedHistogram(xy, nBins);

% plot
figure(1); clf
plotUnityHistogram(xy2, 'nBins', bcenters, 'Offset', offset, 'UnityLine', 0, 'MaxCount', max(count), 'Yscale', 10); hold on
plotUnityHistogram(xy1, 'nBins', bcenters, 'Offset', offset, 'MaxCount', max(count), 'Yscale', 10); hold on
