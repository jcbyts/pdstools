function [plxtrialstart, plxtrialstop] = pdsTrialTimes(PDSname, strobed, events)
% [plxtrialstart, plxtrialstop] = plx_pdsTrialTimes(PDSname, strobed, events)
% get start and stop time for each pds trial in plexon time

if isstruct(PDSname)
    PDS = PDSname;
    clear PDSname;
else
    PDS = [];
    % load PDS file
    fprintf('Loading PDS file\n%s\n', PDSname)
    evalc('load(PDSname, ''-mat'')'); % load PDS file (name entered at top of the file)
end

% find whether pds or plx started first
pdsstart = find(sum(bsxfun(@eq, strobed.values(1,2:end) , PDS.unique_number(:,2:end)),2) == 5);
plxstart = find(sum(bsxfun(@eq, strobed.values(:,2:end) , PDS.unique_number(1,2:end)),2) == 5);
pdsend	 = find(sum(bsxfun(@eq, strobed.values(end,2:end) , PDS.unique_number(:,2:end)),2) == 5);
if isempty(pdsend)
    plxend   = find(sum(bsxfun(@eq, strobed.values(:,2:end) , PDS.unique_number(end,2:end)),2) == 5);
    pdsend = plxend + 1;
else
    plxend   = find(sum(bsxfun(@eq, strobed.values(:,2:end) , PDS.unique_number(pdsend,2:end)),2) == 5);
end

if isempty(pdsstart),
    disp('plx file was started first')
    pdsstart = 1;
else
    disp('pds file was started first')
    plxstart = 1;
end
fprintf('starting match from pds trial %0.0f, plx trial %0.0f \n', pdsstart, plxstart);

ntrials = min(plxend-plxstart,pdsend-pdsstart);

% find trial start event
stoffset = events.time-strobed.times(plxstart);
% stoffset = events.time-strobed.times(3);
% stoffset(stoffset < 0) = inf; % times before strobe don't count
[~, id] = min(abs(stoffset));
trialStartId = events.id(id);
trialStartEventName = events.name{trialStartId};
fprintf('Bit ''%s'' was flipped %d ms from the first strobe. Using that to align trial times\n', trialStartEventName, round(1e3*stoffset(id)))
trialStarts = events.time(events.id==trialStartId);
if ~isfield(PDS, 'timing') || ~isfield(PDS.timing, 'syncTimeDuration') || stoffset(id) > PDS.timing.syncTimeDuration(1) 
    fprintf('something is wrong with your bit flip timing. Using strobe times instead. These are less precise.\n')
    trialStarts = strobed.times;
end


plxtr = plxstart;
pdstr = pdsstart;
plxtrialstart = nan(numel(PDS.goodtrial),1);
plxtrialstop  = nan(numel(PDS.goodtrial),1);

for tr = 1:ntrials
    
    uniqchk = sum(strobed.values(plxtr,2:end) == PDS.unique_number(pdstr,2:end))==5;

    while ~uniqchk
        plxtimestamp = datenum(double([PDS.unique_number(pdstr,1) strobed.values(plxtr,2:end)]));
        pdstimestamp = datenum(PDS.unique_number(pdstr,1:end));
        if plxtimestamp<pdstimestamp
            plxtr = plxtr + 1;
        else
            pdstr = pdstr + 1;
        end
        
        uniqchk = sum(strobed.values(plxtr,2:end) == PDS.unique_number(pdstr,2:end))==5;
    end
    [~, id] = min(abs(strobed.times(plxtr)-trialStarts));
    
    trialstarttime = trialStarts(id);
    if id < numel(trialStarts)
        trialstoptime  = trialStarts(id+1);
    else
        trialstoptime  = max(structfun(@max, plx.events));
    end
   
    plxtrialstart(pdstr) = trialstarttime;
    plxtrialstop(pdstr)  = trialstoptime;
    
    pdstr = pdstr + 1;
    plxtr = plxtr + 1;
    
end
