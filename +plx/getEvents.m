function [events, strobed, info] = plx_getEvents(plxname, info, verbose)
% [events, strobed, info] = plx_getEvents(plxname, info, verbose)
% Get strobe and event timestamps from *.plx file
% inputs:
%   plxname     [string] - full/path/to/*.plx
%   info        [struct] - meta info about the experiment
%       .year - required to parse strobe values
%       
%   Can be run with no inputs. It will open a ui box for selecting the
%   *.plx file. if info is not passed in plx_getInfo will be run to get
%   info.
% 
% output:
%   events  [struct] - fieldnames are the events found in plxname.plx 
%                      each field is a vector of event times.
%   strobed [struct] - 
%        .time  [n x 1] vector of strobe times
%        .value [n x m] matrix of n strobed words
%   info    [struct] - returned from plx_getInfo if it is not passed in
% 

% 20130531  jly     wrote it

if nargin < 3
    verbose = 1;
    if nargin < 2
        info = [];
        if nargin < 1
            [plxfile, plxpath] = uigetfile('*.plx');
            plxname = fullfile(plxpath, plxfile);
        end
    end
end

if isempty(info) || ~isfield(info, 'year')
    info = plx_getInfo(plxname);
end



%-------------------------------------------------------------------------%
% strobed words
[~, timestamps, strobed_values] = plx_event_ts(plxname, 257);
strobe_index   = (strobed_values == mod(info.year,256));
strobe_times   = timestamps(strobe_index);
strobed_values = reshape(strobed_values, 6, [])';

strobed.times  = strobe_times;
strobed.values = strobed_values;

%-------------------------------------------------------------------------%
% find channels that have spike events on them
[~, ~, evcounts] = plx_info(plxname,1);

ch_counts = find(evcounts>1); % event channels

[~,event_names] = plx_event_names(plxname);
event_names = event_names(ch_counts,:);

%-------------------------------------------------------------------------%
% Pull out events
% the plexon file has meta data that defines the start and end of each
% trial as well as the timing of specific trial events
% events.time [n x 1] vector of times
% events.id   [n x 1] vector of index into events.name
% events.name [m x 1] cell array of strings
events = struct();
nEvents =length(ch_counts);
names = cell(1);
if verbose
    fprintf('Pulling out events from the plexon file\r')
end
times = [];
id    = [];
a = 1;
for ii = 1:nEvents
    if ch_counts(ii)<300
        [~, ts] = plx_event_ts(plxname, ch_counts(ii)+1); % weird issue with 0 and 1 base... plexon seems to do this a lot
        % prune Strobed
        strname = (event_names(ii,:));
        whitespace = strfind(strname, ' ');
        fieldname  = strname(setdiff(1:size(strname,2), whitespace));
        
        if ~strcmp(fieldname(1:7), 'Strobed') % strobed values are kept separately
            times = [times; ts(:)];
            id    = [id; a*ones(numel(ts),1)];
            names{a} = fieldname;
            a = a+1;
        end
    end
end

[events.time, order] = sort(times);
events.id = id(order);
events.name = names;



% old format: event was a struct where each field was an event name whos
% values were a vector of times
% %-------------------------------------------------------------------------%
% % Pull out events
% % the plexon file has meta data that defines the start and end of each
% % trial as well as the timing of specific trial events
% events = struct();
% nEvents =length(ch_counts);
% ts = cell(nEvents,1);
% sv = cell(nEvents,1);
% if verbose
%     fprintf('Pulling out events from the plexon file\r')
% end
% 
% for ii = 1:nEvents
%     if ch_counts(ii)<300
%         [~, ts{ii}, sv{ii}] = plx_event_ts(plxname, ch_counts(ii)+1); % weird issue with 0 and 1 base... plexon seems to do this a lot
%         strname = (event_names(ii,:));
%         whitespace = strfind(strname, ' ');
%         fieldname  = strname(setdiff(1:size(strname,2), whitespace));
%         if verbose
%             disp(fieldname)
%         end
%         if ~strcmp(fieldname(1:7), 'Strobed') % strobed values are kept separately
%             events.(fieldname) = ts{ii};
%         end
%     end
% end



