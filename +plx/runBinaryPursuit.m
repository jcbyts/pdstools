function bpspikes = runBinaryPursuit(plxname, spikes, info, greediness, saveDir)
% spikes = plx_runBinaryPursuit(plxname, spikes, info, greediness)

% TODO: set this up to run based on how much memory you have on your
% machine
secondsPerPillowBlock = 10*4*60;%best in multiples of 4 minutes

if ~exist('saveDir', 'var')
    saveDir = pwd;
end

SNRthresh =3;
%% -----------------------------------------------------------------------%
% Binary pursuit: 
% Binary pursuit is an algorithm for detecting simultaneous spikes from
% units on the same channel. You must specify the channel.
% NOTE: this only works for one channel at a time. It
% can be run on tetrode style recordings, but this code will not set
% that up for you. 
if ~isfield(info, 'binaryPursuitChannels') && info.trodalness == 1
    
    ChannelsWithIsolatedUnits = spikes.channel(spikes.snr > SNRthresh);
    
    possibleChannels = unique(ChannelsWithIsolatedUnits);
    info.binaryPursuitChannels = possibleChannels(sum(bsxfun(@eq, possibleChannels, ChannelsWithIsolatedUnits'))>=2);
    if isempty(info.binaryPursuitChannels)
        fprintf('no channels are eligible for binary pursuit algorithm.\nThere must be more that one unit with an spikes.snr > 4 on a single channel to run binary pursuit\n') 
        info = rmfield(info, 'binaryPursuitChannels');
        bpspikes = spikes;
        return
    else
%         yn = input(sprintf('found %d channel(s) eligible for binary pursuit.\nDo you want to run it now (it takes ~1hr per channel)? (y/n)', numel(info.binaryPursuitChannels)), 's');
%         if ~strcmp(yn, 'y')
%             fprintf('skipping binary pursuit\n')
%             info = rmfield(info, 'binaryPursuitChannels');
%         else
            fprintf('running binary pursuit. prepare to wait.\n')
%         end
    end
end
%%


% initialize spike times
final_spike_id      = [];
final_spike_times   = [];
final_spike_waves   = [];
% check if binary pursuit should be run
if isfield(info, 'binaryPursuitChannels') && ~isempty(info.binaryPursuitChannels)
    % loop over channels and run only on specified channels (it takes time)
    for iCh = 1:numel(info.binaryPursuitChannels)
        unitsOnChannel = info.binaryPursuitChannels(iCh) == spikes.channel;
        unitsAreIsolated = spikes.snr > SNRthresh;
        unitsToPursue    = find(unitsOnChannel & unitsAreIsolated);
        % if there are more than one unit on a channel that have an
        % isolated waveform, they can be passed into the binary pursuit
        % algorithm. We have to figure out what to do with multi-unit (MU)
        % activity. It might be that MU is useful for something, but it
        % force this analysis to have to many false positives. I'm not sure
        % about that.
        %%

            spikeIdx = false(size(spikes.time)); 
            for ii = 1:numel(unitsToPursue)
                spikeIdx = spikeIdx | spikes.id==unitsToPursue(ii);
            end
            
            tmp_spike_times = spikes.time(spikeIdx);
            tmp_spike_id    = spikes.id(spikeIdx);
            
            [an_info, an_data] = plx_getAnalog(plxname, info.binaryPursuitChannels(iCh));

            % convert the unit id into a relative number
            C0 = tmp_spike_id;
            unitIds = unique(C0); C = C0;
            for jj  =1:numel(unitIds)
                C(C0==unitIds(jj)) = jj;
            end
            
            % convert spike times from seconds into sample indices
            spike_sample_indices = convertTimeToSamples(tmp_spike_times, info.sampling_rate, an_info.timestamps, an_info.fragsamples);
            
            % reshape into sparse matrix [samples x units]
            spike_matrix = sparse(spike_sample_indices, C, 1, an_info.nsamples, numel(unitIds));
            
            [spikewaves, labels] = pullOutSpikesFromData(an_data,spike_matrix,info.samples_pre_thresh,info.samples_per_waveform,1);
    
            %% plot spike waveforms clipped in Offline Sorter and compare to ones clipped in matlab
            nUnitsToPursue = numel(unitsToPursue);
            wTime = 1e3*((1:info.samples_per_waveform)-info.samples_pre_thresh)/info.sampling_rate;
            figure(33); clf
            cmap = lines;
            for ii = 1:nUnitsToPursue
                % plot clipped waveforms from plx file
                subplot(3,nUnitsToPursue, ii)
                idx = find(spikes.id==unitsToPursue(ii));
                plot(info.waveform_time*1e3, spikes.waveform(idx(1:100:end),:), 'Color', .5*[1 1 1]); hold on
                plot(info.waveform_time*1e3, mean(spikes.waveform(idx, :)), 'Color', cmap(ii,:), 'Linewidth', 2); axis tight
                xlabel('time')
                ylabel('micro volts?')
                title(sprintf('Clipped by Plexon (un: %d, ch: %d, spikes.snr: %02.2f)', unitsToPursue(ii), spikes.channel(unitsToPursue(ii)), spikes.snr(unitsToPursue(ii))))
                % clip your own waveforms with using spike times
                subplot(3, nUnitsToPursue, ii + nUnitsToPursue)
                idx = find(labels==ii);
                plot(wTime, spikewaves(idx(1:100:end),:), 'Color', .5*[1 1 1]); hold on
                plot(wTime, mean(spikewaves(idx,:)), 'Color', cmap(ii,:), 'Linewidth', 2); axis tight
                title(sprintf('Clipped in Matlab (un: %d, ch: %d, spikes.snr: %02.2f)', unitsToPursue(ii), spikes.channel(unitsToPursue(ii)), spikes.snr(unitsToPursue(ii))))
                xlabel('time')
                ylabel('why aren''t the units the same?')
                
            end
            
            %% plot the continuous trace with spike waveforms labeled in color
            subplot(3,nUnitsToPursue,(1:nUnitsToPursue)+nUnitsToPursue*2); 
            hold off
            idx = randi(an_info.nsamples-12e3, 1)+(0:12e3);
            tt = convertSamplesToTime(idx, an_info.adfreq, an_info.timestamps, an_info.fragsamples);
            plot(tt, an_data(idx), 'k'); axis tight; hold on
            spind = spike_sample_indices>idx(1) & spike_sample_indices < idx(end);
            for uu = 1:numel(unitsToPursue)
                spikees = find(spind & C == uu);
                for ii = 1:numel(spikees)
                    jj = spike_sample_indices(spikees(ii))-info.samples_pre_thresh + (1:info.samples_per_waveform);
                    plot(convertSamplesToTime(jj, an_info.adfreq, an_info.timestamps, an_info.fragsamples), an_data(jj), '-', 'Color', cmap(uu,:))
                end
            end
            title('Waveforms labeled on the continuous trace')
            xlabel('time')
            ylabel('what units are these?')
            set(gcf, 'papersize', [8 5], 'paperposition', [0 0 8 5])
            saveas(gcf, fullfile(saveDir, sprintf('binary_pursuit_wf_before_ch%02.0fg_%03.0f.pdf', info.binaryPursuitChannels(iCh), 1/greediness)));
            
            
            %% ------------------------------------------------------------%
            % Run Binary pursuit. This is slow. Don't plan on sitting and
            % waiting for it to run. It takes about an hour for a single
            % channel.
            
            [finalSpikeTimes] = runThePillowAlgorithm(an_data,C,spike_sample_indices, ...
                info.sampling_rate*secondsPerPillowBlock,info.sampling_rate,...
                1e3*info.samples_per_waveform/info.sampling_rate, greediness);

            save(fullfile(saveDir, sprintf('tmp_bpspikes_ch_%03.0f_g_%03.0f.mat', info.binaryPursuitChannels(iCh), 1/greediness)), 'finalSpikeTimes');
            
            %% -----------------------------------------------------------%
            % plot the output
            tmpFinalSpikeTimes = finalSpikeTimes; 
            nNeurons = size(finalSpikeTimes,2);
            % remove short ISI spikes (these can be an artifact of binary
            % pursuit)
            fprintf('Checking for short inter-spike-intervals.\nThese are artifacts of the binary pursuit algorithm and will be removed\n...\n')
            for iNeuron = 1:nNeurons
                spikeSampleTime = find(tmpFinalSpikeTimes(:,iNeuron));
                spikeIsis = diff(spikeSampleTime)/info.sampling_rate/2;
                shortIsis = find(spikeIsis<800e-6); % spikes less than 800 microseconds
                fprintf('neuron %d had %d short Isis\n', iNeuron, length(shortIsis))
                spikesToRemove = spikeSampleTime(shortIsis);
                tmpFinalSpikeTimes(spikesToRemove,iNeuron) = 0;
            end


            % convert back to spike times
            [samples, id] = find(tmpFinalSpikeTimes);
            sptimes = convertSamplesToTime(samples, an_info.adfreq, an_info.timestamps, an_info.fragsamples);
            final_spike_times = [final_spike_times; sptimes(:)];
            final_spike_id    = [final_spike_id; unitsToPursue(id)'];
            
            %%
            [spikees, labels] = pullOutSpikesFromData(an_data,tmpFinalSpikeTimes,info.samples_pre_thresh,info.samples_per_waveform,1);
            final_spike_waves = [final_spike_waves; spikees];
            figure(34); clf
            for ii = 1:nUnitsToPursue
                % clip your own waveforms with using spike times
                subplot(2, nUnitsToPursue, ii)
                idx = find(labels==ii);
                plot(wTime, spikees(idx(1:100:end),:), 'Color', .5*[1 1 1]); hold on
                plot(wTime, mean(spikees(idx,:)), 'Color', cmap(ii,:), 'Linewidth', 2); axis tight
                title(sprintf('(un: %d, ch: %d, snr: %02.2f)', unitsToPursue(ii), spikes.channel(unitsToPursue(ii)), spikes.snr(unitsToPursue(ii))))
                xlabel('time')
                ylabel('units?')
            end
            
            %% plot the continuous trace with spike waveforms labeled in color
            subplot(2,nUnitsToPursue,(1:nUnitsToPursue)+nUnitsToPursue); hold off
            idx = randi(an_info.nsamples-12e3, 1)+(0:12e3);
            tt = convertSamplesToTime(idx, an_info.adfreq, an_info.timestamps, an_info.fragsamples);
            plot(tt, an_data(idx,:), 'k'); axis tight; hold on
            spind = samples>idx(1) & samples < idx(end);
            for uu = 1:numel(unitsToPursue)
                spikees = find(spind & id == uu);
                for ii = 1:numel(spikees)
                    jj = samples(spikees(ii))-info.samples_pre_thresh + (1:info.samples_per_waveform);
%                     plot(time(jj), yy(jj), '-', 'Color', cmap(uu,:))
                    plot(convertSamplesToTime(jj, an_info.adfreq, an_info.timestamps, an_info.fragsamples), an_data(jj), '-', 'Color', cmap(uu,:))
                end
            end
            title('Waveforms labeled on the continuous trace')
            xlabel('time')
            ylabel('what units are these?')
            set(gcf, 'papersize', [8 5], 'paperposition', [0 0 8 5])
            saveas(gcf, fullfile(saveDir, sprintf('binary_pursuit_wf_after_ch%02.0fg_%03.0f.pdf', info.binaryPursuitChannels(iCh), 1/greediness)));
            
            %% -----------------------------------------------------------%
            % plot cross correlation between cells to see the effect of
            % binary pursuit on fine timescale correlation
            cellPairs = nchoosek(1:nUnitsToPursue, 2);
            npairs = size(cellPairs,1);
            figure(5); clf
            for iPair = 1:npairs
                subplot(1, npairs, iPair)
                c1 = cellPairs(iPair,1);
                c2 = cellPairs(iPair,2);
                tb  = 0:1e-3:(size(spike_matrix,1)/info.sampling_rate);
                % get spike time correlation BEFORE binary pursuit
                st1 = find(spike_matrix(:,c1))/info.sampling_rate;
                st2 = find(spike_matrix(:,c2))/info.sampling_rate;
                sb1 = histc(st1, tb);
                sb2 = histc(st2, tb);
                [xc, lags] = xcorr(sb1, sb2, 50, 'coeff');
                plot(lags, xc, 'k'); title(sprintf('xcorr %d, %d (before)', unitsToPursue(c1), unitsToPursue(c2)))
                xlabel('lags(ms)')
                ylabel('corr coef')
                hold on
                % get spike time correlation AFTER binary pursuit
                st1 = find(tmpFinalSpikeTimes(:,c1))/info.sampling_rate;
                st2 = find(tmpFinalSpikeTimes(:,c2))/info.sampling_rate;
                sb1 = histc(st1, tb);
                sb2 = histc(st2, tb);
                [xc, lags] = xcorr(sb1, sb2, 50, 'coeff');
                plot(lags, xc, 'r'); title(sprintf('cross corr %d, %d', unitsToPursue(c1), unitsToPursue(c2)))
                xlabel('lags(ms)')
                ylabel('corr coef')
                legend({'before', 'after'}, 'Location', 'Best')
                set(5, 'papersize', [3*nUnitsToPursue 3], 'paperposition', [0 0 3*npairs 3])
                saveas(5, fullfile(saveDir, sprintf('binary_pursuit_xcorr_ch%02.0fg_%03.0f.pdf', info.binaryPursuitChannels(iCh), 1/greediness)));
            end
            
            %-------------------------------------------------------------%
            % figure 6: autocorrelation for each cell
            figure(6); clf
            for iCell = 1:nUnitsToPursue
                subplot(1, nUnitsToPursue, iCell)
                  st1 = find(spike_matrix(:,iCell))/info.sampling_rate;
                  sb1 = histc(st1, tb);
                  [xc, lags] = xcorr(sb1, sb1, 100, 'coeff');
                  xc(lags==0) = 0;
                  plot(lags, xc, 'k'); title(sprintf('autocorr %d', unitsToPursue(iCell))); hold on
                  st1 = find(tmpFinalSpikeTimes(:,iCell))/info.sampling_rate;
                  sb1 = histc(st1, tb);
                  [xc, lags] = xcorr(sb1, sb1, 100, 'coeff');
                  xc(lags==0) = 0;
                  plot(lags, xc, 'r'); title(sprintf('autocorr %d', unitsToPursue(iCell)))
                xlabel('lags(ms)')
                ylabel('corr coef')
                legend({'before', 'after'}, 'Location', 'Best')
            end
            set(6, 'papersize', [3*nUnitsToPursue 3], 'paperposition', [0 0 3*nUnitsToPursue 3])
            saveas(gcf, fullfile(saveDir, sprintf('binary_pursuit_autocorr_ch%02.0fg_%03.0f.pdf', info.binaryPursuitChannels(iCh), 1/greediness)));
                  
                


    end
end

%% cleanup and build new spikes struct
unitsNoBp = setdiff(unique(spikes.id),unique(final_spike_id));
nUnitsToCombine = numel(unitsNoBp);
for ii = 1:nUnitsToCombine
    idx = spikes.id == unitsNoBp(ii);
    final_spike_id = [final_spike_id; spikes.id(idx)];
    final_spike_times = [final_spike_times; spikes.time(idx)];
    final_spike_waves = [final_spike_waves; spikes.waveform(idx,:)];
end

[sptimes, sortid] = sort(final_spike_times);
bpspikes.time = sptimes;
bpspikes.id   = final_spike_id(sortid);
bpspikes.waveform = final_spike_waves(sortid,:);
bpspikes.channel  = spikes.channel(unique(bpspikes.id));
uniqs = unique(bpspikes.id);
nUnits = numel(uniqs);
for ii = 1:nUnits
    un = uniqs(ii);
    n = sum(bpspikes.id == un);
    spike_waves = bpspikes.waveform(bpspikes.id == un,:);
    mw = mean(spike_waves);
    peak   = max(mw);
    trough = min(mw);
    amp = peak - trough;
    mwbar = repmat(mw, n, 1);
    r = spike_waves - mwbar;
    sigma = std(r(:));
    bpspikes.snr(ii) = amp/(2*sigma);
end


bpspikes.first_continuous_channel = spikes.first_continuous_channel;
bpspikes.continuous_only = spikes.continuous_only;
