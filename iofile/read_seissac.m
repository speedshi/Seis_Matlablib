function seismic=read_seissac(file_seismic,comp)
% This function is used to read seismic data of SAC format.
%
% The seismic data of different traces must have the same sampling rate!
% The different traces may have different starting times and length;
% however in this situation, we will have to cut the traces to let them
% have the same staring time and length.
%
% The input parameter 'file_seismic' can be a string or a cell array which
% contains the names of all the SAC files.
%
% INPUT--------------------------------------------------------------------
% file_seismic: file name of input seismic data;
% comp: character, component of seismic data to read, can be 'Z','N','E';
%
% OUTPUT-------------------------------------------------------------------
% seismic: matlab structure, contains seismic data and the related infomation;
% seismic.network: string, the name of the network;
% seismic.component: character, the name of the data component (usually N, E or Z);
% seismic.name: cell array, 1*ns, contains the name of each station;
% seismic.fe: scaler, the sampling frequency of the data, in Hz;
% seismic.dt: scaler, the time sample interval of the data, in second;
% seismic.t0: matlab datetime, the origin time of the seismic data;
% seismic.NPTS: total number of samples;
% seismic.zdata: 2D array, ns*nt, contains seismic data of Z component;
% seismic.ndata: 2D array, ns*nt, contains seismic data of North component;
% seismic.edata: 2D array, ns*nt, contains seismic data of East component;


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
for ii = 1:n_file
    try
        [~,t0,header]=rdsac(file{ii});
        seismic.dt = header.DELTA; % time sample interval, in second
        seismic.fe = 1.0/header.DELTA; % sampling frequency in Hz (the reciprocal of sampling interval)
        
        if ~exist('time_1','var')
            time_1 = datetime(t0,'ConvertFrom','datenum');  % begin time
            time_2 = time_1 + seconds(header.DELTA*(header.NPTS-1));  % end time
        else
            temp_1 = datetime(t0,'ConvertFrom','datenum');  % begin time
            temp_2 = temp_1 + seconds(header.DELTA*(header.NPTS-1));  % end time
            time_1 = max(time_1,temp_1);
            time_2 = min(time_2,temp_2);
        end
        
    catch ME
        if (strcmp(ME.message,'FILENAME must be a valid file name.'))
            continue;
        else
            rethrow(ME)
        end
    end
end


if ~exist('time_1','var')
    rethrow(ME)
end

seismic.t0 = time_1; % starting time of all the traces
seismic.NPTS = round(seconds(time_2-time_1)*seismic.fe)+1; % total number of time samples

for ii=1:n_file
    % loop through all the files
    try
        [data_temp,t0,header]=rdsac(file{ii});
        seismic.network{ii}=header.KNETWK; % network name of the array
        seismic.component{ii}=header.KCMPNM; % component name of the data
        seismic.name{ii}=header.KSTNM; % name of the station
        temp = datetime(t0,'ConvertFrom','datenum');
        
        id1 = round(seconds(seismic.t0-temp)*seismic.fe)+1;
        id2 = id1 + seismic.NPTS - 1;
        if strcmp(seismic.component{ii}(end),'Z') || strcmp(seismic.component{ii}(end),'z')
            seismic.zdata(ii,:) = data_temp(id1:id2);
        elseif strcmp(seismic.component{ii}(end),'N') || strcmp(seismic.component{ii}(end),'n')
            seismic.ndata(ii,:) = data_temp(id1:id2);
        elseif strcmp(seismic.component{ii}(end),'E') || strcmp(seismic.component{ii}(end),'e')
            seismic.edata(ii,:) = data_temp(id1:id2);
        end
    catch ME
        if (strcmp(ME.message,'FILENAME must be a valid file name.'))
            % no data
            seismic.network{ii} = 'NULL';
            seismic.component{ii} = 'NULL';
            seismic.name{ii} = 'NULL';
            if strcmp(comp(end),'Z') || strcmp(comp(end),'z')
                seismic.zdata(ii,:) = NaN(seismic.NPTS,1);
            elseif strcmp(comp(end),'N') || strcmp(comp(end),'n')
                seismic.ndata(ii,:) = NaN(seismic.NPTS,1);
            elseif strcmp(comp(end),'E') || strcmp(comp(end),'e')
                seismic.edata(ii,:) = NaN(seismic.NPTS,1);
            end
            
            continue;
        else
            rethrow(ME)
        end
    end
end


end