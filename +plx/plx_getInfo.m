function info = plx_getInfo(plxname, info)
% info = plx_getInfo(plxname, info)
% Get basic info about plexon file
% inputs:
%   plxname     [string] - full/path/to/*.plx
%   info        [struct] - meta info about the experiment
%   Can be run with no inputs. It will open a ui box for selecting the
%   *.plx file.
% 
% output:
%   info is updated to contain the fields
%       .full_path_to_file      - the name of the plexon file
%       .sort_client_version    - version of sort client
%       .samplesPerWaveform     - number of samples per spike waveform
%       .sampling_rate          - sampling rate for clipped waveforms
%       .samplesPreThresh       - number of samples before threshold crossing
%       .comment
%       .trodalness
%       .file_duration
%       .samples_pre_thresh
%       .samples_per_waveform


if nargin < 2
    info = [];
    
    if nargin < 1
        [plxfile, plxpath] = uigetfile('*.plx');
        plxname = fullfile(plxpath, plxfile);
    end
end


[info.full_path_to_file, info.sort_client_version, info.sampling_rate, info.comment, info.trodalness, info.samples_per_waveform, info.samples_pre_thresh, info.spike_peak_mV, info.spike_res_bits, info.analog_peak_mV, info.analog_res_bits, info.file_duration, DateTime] = plx_information(plxname);

% parse date-time. 'year' is used to extract unique trial words
sid = strfind(DateTime, '/');
info.month = str2double(DateTime(1:sid(1)-1));
info.day   = str2double(DateTime(sid(1)+1:sid(2)-1));
info.year  = str2double(DateTime(sid(2)+(1:4)));
cid = strfind(DateTime, ':');
info.hour = str2double(DateTime(cid(1)-2+(0:1)));
info.minutes = str2double(DateTime(cid(1)+(1:2)));
info.seconds = str2double(DateTime(cid(2)+(1:2)));

info.waveform_time = (-info.samples_pre_thresh:(info.samples_per_waveform-info.samples_pre_thresh-1))/info.sampling_rate;