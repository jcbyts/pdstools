function [events, strobed] = getEvents(pl, info, verbose)
% Get strobe and event timestamps from *.plx file
% [events, strobed, info] = getEvents(plstruct, info, verbose)
% inputs:
%    plstruct   [struct] - output of readPLX
%       
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
% 20140820  jly     replaced with readPLX version

import plx.*

if nargin < 2
    verbose = 1;
    if nargin < 1
        help plx.getEvents
        return
    end
end

%-----------------------------------------------------------------------------%
% make strobed struct
V = datevec(pl.Date);
%
eventNames = {pl.EventChannels(:).Name};
strobeChannel = cellfun(@(x) strcmp(x, 'Strobed'), eventNames);
assert(strcmp(pl.EventChannels(strobeChannel).Name, 'Strobed'))
pl.Year = V(1);

strobeIndex = pl.EventChannels(strobeChannel).Values == mod(pl.Year, 256);
strobed.times  = double(pl.EventChannels(strobeChannel).Timestamps(strobeIndex))/pl.ADFrequency;
strobed.values = reshape(pl.EventChannels(strobeChannel).Values, 6, [])';

%-----------------------------------------------------------------------------%
% make events struct
nEvents = numel(pl.EventChannels)-1;
events.time = [];
events.id   = [];
events.name = cell(1);
eventCtr = 1;
for ii = 1:nEvents
    if pl.EventChannels(ii).Channel > 1 && pl.EventChannels(ii).Channel < 256
        events.time = [events.time; double(pl.EventChannels(ii).Timestamps)/pl.ADFrequency];
        events.id   = [events.id; ones(numel(pl.EventChannels(ii).Timestamps),1)*eventCtr];
        events.name{eventCtr} = pl.EventChannels(ii).Name;
        eventCtr = eventCtr + 1;
    end
end

[events.time, ord] = sort(events.time);
events.id = events.id(ord);
