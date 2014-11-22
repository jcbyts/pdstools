function pdsCombine(plxname, pdsname, info, verbose, useContinuousOnly)
% pdsCombine(plxname, pdsname, info, verbose, useContinuousOnly)
% DEPRECIATED??

if nargin < 5
    useContinuousOnly = 0;
    if nargin < 4
        verbose = 1;
        if nargin < 3
            info = [];
            if nargin < 2
                [pdsfile, pdspath] = uigetfile('*.PDS');
                pdsname = fullfile(pdspath, pdsfile);
                if nargin < 1
                    [plxfile, plxpath] = uigetfile('*.plx');
                    plxname = fullfile(plxpath, plxfile);
                end
            end
        end
    end
end

[saveDir, fname] = fileparts(pdsname, 'all');

pl = readPLXFileC(plxname);

V = datevec(pl.Date);
assert(strcmp(pl.EventChannels(end).Name, 'Strobed'))
year = V(1);

strobe_index = pl.EventChannels(end).Values == mod(year, 256);

%-------------------------------------------------------------------------%
% step 1: get general info
info = plx_getInfo(plxname, info);

%-------------------------------------------------------------------------%
% step 2: get event data
[events, strobed, info] = plx_getEvents(plxname, info, verbose);

%-------------------------------------------------------------------------%
% step 3: get spike data from offline sorter
spikes = plx_getSpikes(plxname, useContinuousOnly, info, verbose);

%% ---------------------------------------------------------------------------%

savename = fullfile(saveDir, [fname '_events.mat']);
if ~exist(savename, 'file')
    save(savename, 'strobed', 'events', '-v7.3')
end
savename = fullfile(saveDir, [fname '_info.mat']);
if ~exist(savename, 'file')
    save(savename, 'info', '-v7.3')
end

% get pds trial starts
[start, stop] = plx_pdsTrialTimes(pdsname, strobed, events);

% save(fullfile(saveDir, [fname '_events.mat']), 'strobed', 'events', 'info', '-v7.3')

spikes = assignTrialSpikes(spikes, start, stop);

%-------------------------------------------------------------------------%
% step 4: get LFP
% find LFP channels
savename = fullfile(saveDir, [fname '_lfp.mat']);
if exist(savename, 'file')
    fprintf('LFP file already exists. Skipping\n\n\n\n')
else
    [~, adfreq] = plx_adchan_freqs(plxname);
    [~, adc]    = plx_ad_chanmap(plxname);
    lfp_samplingrate = min(adfreq);
    lfp_channels     = adc(adfreq==lfp_samplingrate)+1; % + 1 because we like base 1
    
    % pull out data 
    [lfp_info, lfp_data] = plx_getAnalog(plxname, lfp_channels);
    
    lfp_info.trial_start = convertTimeToSamples(start, lfp_info.adfreq, lfp_info.timestamps, lfp_info.fragsamples);
    lfp_info.trial_stop  = convertTimeToSamples(stop,  lfp_info.adfreq, lfp_info.timestamps, lfp_info.fragsamples);
    fprintf('saving LFP data\n')
    save(savename, 'lfp_info', 'lfp_data', '-v7.3')
end

%% -------------------------------------------------------------------------%
% prepare for Binary pursuit
savename = fullfile(saveDir, [fname '_spikes.mat']);
save(savename, 'spikes', '-v7.3')
if useContinuousOnly
    spikes0 = spikes;
    greeds = 1/1000; %1./[1 10 100 1000];
    for gg = 1:numel(greeds)
        grd = greeds(gg);
        savename = fullfile(saveDir, [fname sprintf('_bp%04.0fspikes.mat', 1/grd)]);
        if exist(savename, 'file')
            fprintf('Binary Pursuit with greed %d already exists. Skipping...\n\n\n\n', 1/grd)
        else
            spikes = plx_runBinaryPursuit(plxname, spikes0, info, grd, saveDir);
            save(savename, 'spikes', '-v7.3')
        end
    end
end
%% -------------------------------------------------------------------------%
% step 6: get analog
% find Raw channels
savename = fullfile(saveDir, [fname '_analog.mat']);
if exist(savename, 'file')
    fprintf('Raw file already exists. Skipping\n')
else
    [~, adfreq] = plx_adchan_freqs(plxname);
    [~, adc]    = plx_ad_chanmap(plxname);
    an_samplingrate = adfreq == info.sampling_rate;
    an_channels     = adc(an_samplingrate)+1; % + 1 because we like base 1
    
    % pull out data
    [an_info, an_data] = plx_getAnalog(plxname, an_channels);
    
    an_info.trial_start = convertTimeToSamples(start, an_info.adfreq, an_info.timestamps, an_info.fragsamples);
    an_info.trial_stop  = convertTimeToSamples(stop,  an_info.adfreq, an_info.timestamps, an_info.fragsamples);
    
    fprintf('saving raw analog data\n')
    save(savename, 'an_info', 'an_data', '-v7.3')
end







