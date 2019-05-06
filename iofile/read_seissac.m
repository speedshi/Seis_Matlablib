function seismic=read_seissac(file_seismic)
% This function is used to read seismic data of SAC format.
%
% The seismic data of different traces should have the same length, the
% same sampling rate and also start at the same time.
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

for ii=1:n_file
    % loop through all the files
    [seismic.data(ii,:),t0,header]=rdsac(file{ii});
    seismic.network{ii}=header.KNETWK; % network name of the array
    seismic.component{ii}=header.KCMPNM; % component name of the data
    seismic.name{ii}=header.KSTNM; % name of the station
    
    if ii==1
       seismic.fe=1./header.DELTA; % sampling frequency in Hz (the reciprocal of sampling interval)
       seismic.t0=datetime(t0,'ConvertFrom','datenum'); % starting time of all the traces
    end
    
end


end