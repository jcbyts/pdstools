function plotUnityHistogram(xy, varargin)
% plot the difference histogram on the unity line
% haven't figured out which options to include yet
% plotUnityHistogram(xy)

import pdsa.*

h=plot(xy(:,1), xy(:,2), 'o'); hold on
axis equal

clr=get(h, 'Color');


v=[1 -1]; v=v/norm(v); % line orthogonal to unity

[bx,by,bcenters, count] = projectedHistogram(xy, 25, v);

% plot histogram
plot(xlim, xlim, 'k') % plot unity line
fill(bx, by, 'k', 'FaceColor', clr)


% plot axes
% plot(bcenters([1 end])*v(1)+off, bcenters([1 end])*v(2)+off, '-')

% plot projected points
% plot(x(:,1)+off, x(:,2)+off, 'or')


% plot(bx, by, '-')
