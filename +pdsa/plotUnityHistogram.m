function h=plotUnityHistogram(xy, varargin)
% plot the difference histogram on the unity line
% haven't figured out which options to include yet
% plotUnityHistogram(xy)
% optional inputs
%   'Color'         [1 x 3] color of plot (default: next color on current axis)
%   'FaceColor'     [1 x 3] color of histogram (default: color of plot)
%   'nBins'         [1 x 1] default 25. TODO: make vector for binedges
%   'Marker'        [char]  (default: 'o')
%   'MarkerSize'    [1 x 1] (default: 10)
%   'MarkerFaceColor' [1 x 3] color of scatter center (default: color of plot)
%   'Linewidth'     [1 x 1] default: 1
%   'UnityLine'     [boolean] plot the unity line (default: true)
%   'Axes'          [char] histogram axes: 'on' or 'off' (default: 'on')
%   'AxesLabels'    [char] histogram axes labels: 'on' or 'off' (default: 'on')
%
% see also: projectedHistogram.m

import pdsa.*

jj = plot(nan,nan);
nextcolor = get(jj,'color');

p=inputParser();
p.addOptional('Color', nextcolor)
p.addOptional('FaceColor', [])
p.addOptional('nBins', 25)
p.addOptional('Marker', 'o')
p.addOptional('MarkerSize', 10)
p.addOptional('MarkerFaceColor', nextcolor);
p.addOptional('Linewidth', 1)
p.addOptional('UnityLine', true)
p.addOptional('Axes', 'on')
p.addOptional('AxesLabels', 'on')
p.parse(varargin{:})

if isempty(p.Results.FaceColor)
    FaceColor=p.Results.Color;
else
    FaceColor=p.Results.FaceColor;
end


v=[1 -1]; v=v/norm(v); % line orthogonal to unity

if strcmpi(p.Results.Axes, 'on')
    [bx,by,bcenters, count] = projectedHistogram(xy, p.Results.nBins, v);
else
    [bx,by]=projectedHistogram(xy, p.Results.nBins, v);
end

xd(1)=min(min(xy));
m=max(max(bx),max(by));
xd(2)=m+.5*(m-min(by));
if p.Results.UnityLine
    plot(xd, xd, 'k') % plot unity line
    hold on
end

h=plot(xy(:,1), xy(:,2), 'o', 'Color', p.Results.Color, 'Marker', p.Results.Marker, 'MarkerFaceColor', p.Results.MarkerFaceColor); hold on




% plot histogram
fill(bx, by, 'k', 'EdgeColor', p.Results.Color, 'FaceColor', FaceColor, 'Linewidth', 1)

% plot axes
if strcmp(p.Results.Axes, 'on')
%     v=v/norm(v);
    u=[1 -1]./v; u=u/norm(u);
    off=([bx(1) by(1)]*u')*u;
    
    % x axis
    bcenters=bcenters(1:2:end);
    inds=[1 numel(bcenters)];

    %     inds=1:numel(bcenters);
    xax=bcenters(inds)*v(1)+off(1);
    xay=bcenters(inds)*v(2)+off(2);
    plot(xax, xay, 'k-', 'Linewidth', 1.5)
    % x ticks
    
    tkx=xax(:)';
    tky=xay(:)';
    tkx=[tkx; tkx-.1*u(1)];
    tky=[tky; tky-.1*u(2)];
    plot(tkx, tky, 'k', 'Linewidth', 1.5)
    
    
    X=([bx by]*u')*u;
    x=[bx(1) by(1)];
    xu=x-(x*u')*u; % orthogonal to u
    X=bsxfun(@plus, X, xu);
    mxX=max(X);
    mnX=min(X);
    Y=[mxX; mnX];
    yax=Y(:,1);
    yay=Y(:,2);
    plot(yax, yay, 'k', 'Linewidth', 1.5)
    
    tkx=yax(:)';
    tky=yay(:)';
    tkx=[tkx; tkx-.1*v(1)];
    tky=[tky; tky-.1*v(2)];
    plot(tkx, tky, 'k', 'Linewidth', 1.5)
    th=cart2pol(v(1), v(2));
    th=th/pi*180;
    for b=1:numel(inds)
        jj=text(xax(b)-.25*u(1),xay(b)-.25*u(2),num2str(bcenters(inds(b)), '%01.1f'));
        set(jj, 'rotation', th)
    end
    
    rcnt=[max(count) min(count)];
    for c=1:2
        jj=text(yax(c)-.5*v(1), yay(c)-.5*v(2), num2str(rcnt(c), '%01.1f'));
        set(jj, 'rotation', th)
    end
end

xlim(xd)
ylim(xd)
axis equal