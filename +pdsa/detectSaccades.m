function [result, smoothTrace] = detectSaccades(sampleTimes,eyeData, varargin)
% [result, smoothTrace] = saccadeDetector(time,trace,varargin)
% Detect saccades from an X,Y trace using velocity thresholds
% Inputs:
%   sampleTimes  [1 x n]    Time of each sample in eyeData
%   eyeData      [2 x n]    XY position of the eye in degrees
% Optional Arguments as argument-pair input (default)
%   'filterPosition' (1)    Smooth the position traces by <n> samples
%   'filterLength'  (20)    Length (samples) of boxcar for smoothing velocity
%   'detectThresh' (100)    Threshold (deg/sec) to detect saccades
%   'startThresh'    (5)    Threshold (deg/sec) for the start of a saccade
%   'minIsi'        (25)    Minimum ISI (samples)
%   'blinkIsi'      (50)    Padding around blinks (samples)
%   'velocityMode' (true)   Project onto saccade vector to measure velocity
%   'verbose',    (false)   Report some details about the saccade detection
% Output:
%   result (struct)
%       .start         saccade start time
%       .end           saccade end time
%       .duration      saccade duration
%       .size          saccade size
%       .startXpos     x position at start
%       .startYpos     y position at start
%       .endXpos       x position at end
%       .endYpos       y position at end
%       .startIndex    start index (into data trace)
%       .endIndex      end index
%   smoothTrace [2 x n]     Smoothed position trace

if size(sampleTimes,1) > size(sampleTimes,2)
    sampleTimes = sampleTimes';
end

ip=inputParser();
ip.addOptional('filterPosition', 1)
ip.addOptional('filterLength', 20);
ip.addOptional('detectThresh', 200);
ip.addOptional('startThresh', 5)
ip.addOptional('minIsi', 25)
ip.addOptional('minDur', 5)
ip.addOptional('blinkIsi', 50)
ip.addOptional('velocityMode', true)
ip.addOptional('verbose', true)
ip.parse(varargin{:});

% --- filter position
pos_filter_length = ip.Results.filterPosition; % short filter to smooth for noisy eye traces

if pos_filter_length > 1
    eyeData = filtfilt(ones(pos_filter_length,1)/pos_filter_length, 1, eyeData')';
end

filter_length = ip.Results.filterLength; % number of amples to average for baseline velocity
detect_thresh = ip.Results.detectThresh; % threshold in deg/s to detect a saccade
start_thresh  = ip.Results.startThresh;  % threshold in deg/s to determine the start and end of a saccade
min_isi       = ip.Results.minIsi;       % minimum number of samples between any two saccades
minDur        = ip.Results.minDur;       % minimum duration of a saccade (samples)
blink_isi     = ip.Results.blinkIsi;     % ignore saccades this many samples around a blink
reproject_onto_saccade_direction   = ip.Results.velocityMode; % work on signed velocity instead of speed

%% calculate velocities
Fs  = mode(diff(sampleTimes)); % find sampling rate
vel = [[0;0],(eyeData(:,3:end)-eyeData(:,1:end-2))/(Fs*2),[0;0]];
speed = sqrt(sum(vel.^2));


if filter_length > 1
    % boxcar filter
    a = 1;
    b = ones(1,filter_length/2)/filter_length/2;
%     speedf = filter(b,a,speed);
    speedf = filtfilt(b,a,speed);
    speedspeedf=speed-speedf;
else
    speedf = speed;
    speedspeedf = speed;
end



% find eye speeds that exceed the detection threshold
potential_saccades = diff([0 (speedspeedf) > detect_thresh]) == 1;
potential_saccades = find(potential_saccades);
% remove isi violations
if numel(potential_saccades) > 1
    potential_saccades([min_isi diff(potential_saccades)]<min_isi)=[];
end

% --- Loop over potential saccades. Estimate some parameters of the saccade
numPotential = length(potential_saccades);
startIndex = zeros(numPotential,1);
endIndex   = zeros(numPotential,1);
baseSpeed  = zeros(numPotential,1);
endSpeed = zeros(numPotential,1);

if ip.Results.verbose
    fprintf('Found %d potential saccades\n', numPotential)
end

if ip.Results.verbose
    figure;
end

for iSaccade = 1:numPotential
%     if ip.Results.verbose
%         iix = -100:500;
%         inds = potential_saccades(iSaccade)+iix;
%         valid = inds > 0 & inds < size(eyeData,2);
%         
%         iix = iix(valid);
%         inds = inds(valid);
%         
%         this_trace = eyeData(:,inds);
%         this_speed = speed(inds);
%         this_speedspeedf = speedspeedf(inds);
%         this_velocity = vel(:,inds);
%         
%         subplot(2,1,1)
%         hold off;
%         plot(iix,this_trace');
%         axis tight
%         
%         hold on;
%         plot([0 0], [min(min(this_trace)) max(max(this_trace))], '-k');
% %         plot([endIndex(iSaccade) endIndex(iSaccade)]-startIndex(iSaccade), [min(min(this_trace)) max(max(this_trace))], '--k');
%         
%         subplot(2,1,2)
%         hold off;
%         plot(iix,this_speed);
%         hold all;
%         plot(iix,this_velocity');
%         plot(iix,this_speedspeedf');
%         
% %         plot([0 0], [min(this_speed) max(this_speed)], '-k');
% %         plot([endIndex(iSaccade) endIndex(iSaccade)]-startIndex(iSaccade), [min(this_speed) max(this_speed)], '--k');
%         
%         axis tight
%         drawnow
%         waitforbuttonpress
%     end
    
    % if this potential saccade occured during the previous saccade. Skip
    if (iSaccade > 1) && (potential_saccades(iSaccade) < endIndex(iSaccade-1))
        startIndex(iSaccade) = potential_saccades(iSaccade);
        endIndex(iSaccade)   = potential_saccades(iSaccade);
        baseSpeed(iSaccade)  = speedf(startIndex(iSaccade));
        continue;
    end
    
    % --- compute parameters
    % 1. Find the first bin where the unfiltered speed minus the filtered 
    %    speed exceeds the threshold
    startIndex(iSaccade) = potential_saccades(iSaccade);
    while speedspeedf(startIndex(iSaccade)) > start_thresh
        startIndex(iSaccade) = startIndex(iSaccade) - 1;
    end
    
    % 2. determine pursuit velocity using the filtered speed trace
    baseSpeed(iSaccade) = speedf(startIndex(iSaccade));
    
    % initialize the end index
    endIndex(iSaccade)  = potential_saccades(iSaccade);
    
    % walk forward sample by sample until the raw speed return below the
    % threshold
    while speed(endIndex(iSaccade)) - baseSpeed(iSaccade) > start_thresh
        endIndex(iSaccade) = endIndex(iSaccade)+1;
    end
    
    if endIndex(iSaccade)+blink_isi > size(speed,2) ... % end exceeds data duration
            || isnan(sum(speed(endIndex(iSaccade) + (1:blink_isi)))) ... % not valid data
            || startIndex(iSaccade) <= blink_isi ... 
            || isnan(sum(speed(startIndex(iSaccade)-(1:blink_isi))))
        % store the start and end index as the same value
        startIndex(iSaccade) = potential_saccades(iSaccade);
        endIndex(iSaccade)   = potential_saccades(iSaccade);
        continue;
    end
    
	if reproject_onto_saccade_direction
    %3.now with these estimators of the actual saccade, we will rebase the
    % velocity along the saccade trajectory to allow using a velocity
    % instead of speed
        
        % get the saccade vector (from start to end point)
        saccade_vector = eyeData(:,endIndex(iSaccade))-eyeData(:,startIndex(iSaccade));
        
        % make unit vector (can rotate only. no scaling)
        expected_basis = saccade_vector/norm(saccade_vector); %sqrt(sum(saccade_vector.^2));
        
        % project on the saccade vector
        vxy=expected_basis'*vel(:,(startIndex(iSaccade)-filter_length):(endIndex(iSaccade)+filter_length));
        
        % recalculate baseline speed along the saccade vector
        baseSpeed(iSaccade) = sum(vxy(1:filter_length))/filter_length;
        
        endSpeed(iSaccade)= sum(vxy(end-filter_length+1:end))/filter_length;
        
        % initialize the old start index
        old_startindex=startIndex(iSaccade);
        
        % find the actual start
        startIndex(iSaccade) = potential_saccades(iSaccade);

        while vxy(1,startIndex(iSaccade)-old_startindex+filter_length+1) - baseSpeed(iSaccade)> start_thresh
            startIndex(iSaccade) = startIndex(iSaccade)-1;
        end
        
        % find the actual end of saccade
        endIndex(iSaccade)=potential_saccades(iSaccade);
        
        % find the actual end
        while vxy(1,endIndex(iSaccade)-old_startindex+filter_length+1) - endSpeed(iSaccade)> start_thresh
            endIndex(iSaccade) = endIndex(iSaccade)+1;
        end
	end
     
end

% if there were no saccades, exit the function
if isempty(iSaccade)
    result = [];
    smoothTrace = [];
    return
end

% --- compute some basic states
sacdur  = endIndex-startIndex;
sacsize = sqrt(sum((eyeData(:,endIndex)-eyeData(:,startIndex)).^2));
%do a hist of saccade durarion, hopefullt a clear cutoff
%  hist(sacdur, 0:350)
%  hist(sacsize,0:60)

% remove saccades below minimum duration
removeIx = sacdur < minDur; 
if ip.Results.verbose
    fprintf('%d potential saccades were shorter than the minimum allowed duration (%d)\n', sum(removeIx), minDur);
end

startIndex(removeIx) = [];
endIndex(removeIx) = [];
baseSpeed(removeIx) = [];
potential_saccades(removeIx) = [];
 
% recompute duration and size
sacdur  = endIndex-startIndex;
sacsize = sqrt(sum((eyeData(:,endIndex)-eyeData(:,startIndex)).^2));

result = struct('start', sampleTimes(startIndex)', ...
    'end', sampleTimes(endIndex)', ...
    'duration', sacdur, ...
    'size', sacsize', ...
    'startXpos', eyeData(1,startIndex)', ...
    'startYpos', eyeData(2,startIndex)', ...
    'endXpos', eyeData(1,endIndex)', ...
    'endYpos', eyeData(2,endIndex)', ...
    'startIndex', startIndex, ...
    'endIndex', endIndex);
    
% result=[sampleTimes(startIndex); sampleTimes(endIndex); sacdur; sacsize; eyeData(:,startIndex); eyeData(:,endIndex); startIndex; endIndex];

% % loop over saccades and plot the analysis
% if ip.Results.verbose
%     figure;
%     for iSaccade= 1:length(potential_saccades)
%         this_trace=eyeData(:,startIndex(iSaccade)+(-100 : 2000));
%         this_speed=speed(startIndex(iSaccade)+(-100 : 2000));
%         this_speedspeedf=speedspeedf(startIndex(iSaccade)+(-100 : 2000));
%         this_velocity = vel(:,startIndex(iSaccade)+(-100 : 2000));
%         %     this_pupil = pupil(startindex(iSaccade)+(-100 : 2000));
%         subplot(2,1,1)
%         hold off;
%         plot(-100:2000,this_trace');
%         hold on;
%         plot([0 0], [min(min(this_trace)) max(max(this_trace))], '-k');
%         plot([endIndex(iSaccade) endIndex(iSaccade)]-startIndex(iSaccade), [min(min(this_trace)) max(max(this_trace))], '--k');
%         subplot(2,1,2)
%         hold off;
%         plot(-100:2000,this_speed);
%         hold all;
%         plot(-100:2000,this_velocity');
%         plot([0 0], [min(this_speed) max(this_speed)], '-k');
%         plot([endIndex(iSaccade) endIndex(iSaccade)]-startIndex(iSaccade), [min(this_speed) max(this_speed)], '--k');
%         
%         waitforbuttonpress
%     end
% end

if nargout > 1
    smoothTrace = eyeData;
end
