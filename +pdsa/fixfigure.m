function [] = fixfigure(figurehandle, fontSize, paperSize, varargin)


if nargin < 3
    paperSize = [6 3.5];
    if nargin < 2
        fontSize = 12;
    end
end


% options (not implemented yet)
p=inputParser();
p.addOptional('TickDir', 'out');
p.addOptional('FontName', 'Helvetica');
p.addOptional('FontWeight','normal');
p.addOptional('LineWidth',.5);
p.addOptional('FontSize', fontSize);
p.parse(varargin{:})
opts=p.Results;
axisOptions= fieldnames(opts);

% Check to see if a specific figure is identified for
% sciencestagramification.  Otherwise the current figure is used:
if exist('figurehandle','var')
    set(0,'currentfigure',figurehandle);
end

figureChildren = get(gcf,'children');

% Loop through all the children of the figure:
for kChild = 1:length(figureChildren)
    currentChildProperties = get(figureChildren(kChild));
    if isfield(currentChildProperties, 'ActivePositionProperty')
        set(gcf,'currentaxes',figureChildren(kChild))
        set(gca, 'YColor', [0 0 0], 'XColor', [0 0 0], 'ZColor', [0 0 0], 'Layer', 'top')
        axisChildren = get(gca,'children');
        currentAxisProperties = get(gca);
        for ii = 1:numel(axisOptions)
%                 disp(axisOptions{ii})
                if isfield(currentAxisProperties, axisOptions{ii}) && strcmpi(currentAxisProperties.Type, 'axes')
                    set(figureChildren(kChild), axisOptions{ii}, opts.(axisOptions{ii}));
                end
        end
        
        % And loop through all the children of each axis:
        for kAxis = 1:length(axisChildren)
            axisfields = get(axisChildren(kAxis));
            
            % Loop through axis options and modify the axis
            for ii = 1:numel(axisOptions)
%                 disp(axisOptions{ii})
                if isfield(axisfields, axisOptions{ii}) && ~strcmpi(axisOptions{ii}, 'Linewidth')
                    
                    set(axisChildren(kAxis), axisOptions{ii}, opts.(axisOptions{ii}));
                end
            end
            
        end
        


        


        ht = get(gca,'title');
        set(ht,'FontName',opts.FontName,'FontSize',opts.FontSize, 'Color', 'k');
        hx = get(gca,'xlabel');
        set(hx,'FontName',opts.FontName, 'FontWeight', opts.FontWeight,'FontSize',opts.FontSize, 'Color', 'k');
        hy = get(gca,'ylabel');
        set(hy,'FontName',opts.FontName, 'FontWeight', opts.FontWeight,'FontSize',opts.FontSize, 'Color', 'k');
        hl  = findobj(gcf,'Type','axes','Tag','legend');
        set(hl,'box','off')
    end
    
    %         % Movshonize the axes by covering the corner
%         xd = xlim;
%         yd = ylim;
%         hold on
%         plot(xd(1), yd(1), '+b', 'markersize', 25, 'Linewidth', 5)
        box off
    
end


set(gcf, 'Papersize', paperSize, 'paperposition', [0 0 paperSize])
set(gcf, 'Color', 'w')


% % This function is used to bring the original axes
% % on the foreground and the transformed axis to the background
% function switch_objects_depth(parenthandle,obj1,obj2);
% children = get(parenthandle,'Children');
% obj1pos = find(children==obj1);
% obj2pos = find(children==obj2);
% children(obj1pos) = obj2;
% children(obj2pos) = obj1;
% set(parenthandle,'Children',children);