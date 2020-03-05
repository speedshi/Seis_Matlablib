function seismic=read_seissac(file_seismic)
% This function is used to read seismic data of SAC format.
%
% The seismic data of different traces must have the same sampling rate.
% The different traces may have different starting times and length;
% however in this situation, we will have to cut the traces to let them
% have the same staring time and length.
%
% The input parameter 'file_seismic' can be a string or a cell array which
% contains the names of all the SAC files.
%
% INPUT--------------------------------------------------------------------
% seismic.network: string, the name of the network;
% file_seismic: file name of input seismic data;
%
% OUTPUT-------------------------------------------------------------------
% seismic: matlab structure, contains seismic data and the related infomation;
% seismic.network: string, the name of the network;
% seismic.component: character, the name of the data component (usually N, E or Z);
% seismic.name: cell array, 1*ns, contains the name of each station;
% seismic.fe: scaler, the sampling frequency of the data, in Hz;
% seismic.dt: scaler, the time sample interval of the data, in second;
% seismic.t0: matlab datetime, the origin time of the seismic data;
% seismic.data: 2D array, ns*nt, contains seismic data.

if isa(file_seismic,'cell')
    % input is cell array which might contain names of different SAC files;
    file=file_seismic;
else
    % input is characters or string which is the name of a SAC file;
    file={file_seismic};
end

n_file=length(file); % the number of SAC files

% check and obtain time and sampling information, keep consistent over all
% traces
[~,t0,header]=rdsac(file{1});
seismic.fe=1.0/header.DELTA; % sampling frequency in Hz (the reciprocal of sampling interval)
seismic.dt=header.DELTA; % time sample interval, in second
time_1 = datetime(t0,'ConvertFrom','datenum');  % begin time
time_2 = time_1 + seconds(header.DELTA*(header.NPTS-1));  % end time
for ii = 2:n_file
    [~,t0,header]=rdsac(file{ii});
    if seismic.dt ~= header.DELTA
        % the sampling interval should be the same
        error('The sampling internal of different traces is different! Check input data!');
    end
    temp_1 = datetime(t0,'ConvertFrom','datenum');
    temp_2 = temp_1 + seconds(header.DELTA*(header.NPTS-1));
    
    time_1 = max(time_1,temp_1);
    time_2 = min(time_2,temp_2);
    
end


seismic.t0 = time_1; % starting time of all the traces
nt = round(seconds(time_2-time_1)*seismic.fe)+1; % number of time samples for each trace


for ii=1:n_file
    % loop through all the files
    [data_temp,t0,header]=rdsac(file{ii});
    seismic.network{ii}=header.KNETWK; % network name of the array
    seismic.component{ii}=header.KCMPNM; % component name of the data
    seismic.name{ii}=header.KSTNM; % name of the station
    temp = datetime(t0,'ConvertFrom','datenum');
    
    id1 = round(seconds(seismic.t0-temp)*seismic.fe)+1;
    id2 = id1 + nt - 1;
    seismic.data(ii,:) = data_temp(id1:id2);
    
end


end