function seismic=read_seish5(fname)
% This function is used to read seismic data of h5 format.
% The length of seismic data of different stations should be the same.
% At the present, it is only valid for reading one component dataset.
% The H5 data should be oganized as follows:
% /NETWORK_NAME (contains different stations)
% --/NETWORK_NAME/STATION_NAME (contains single component data)
% ----/NETWORK_NAME/STATION_NAME/COMPONENT_NAME (contains seismic data)
% /_metadata (contains sampling frequency and time information)
% --/_metadata/fe (sampling frequency)
% --/_metadata/t0_UNIX_timestamp (origin time of the data)
% INPUT--------------------------------------------------
% fname: path and file name of the input seismic data;
% OUTPUT-------------------------------------------------
% seismic: structure, contains seismic data and metadata of each station;
% seismic.network: string, the name of the network;
% seismic.component: character, the name of the data component (usually N, E or Z);
% seismic.name: cell array, contains the name of each station;
% seismic.fe: scaler, the sampling frequency of the data;
% seismic.t0: datetime, the origin time of the seismic data;
% seismic.data: 2D array, ns*nt, contains seismic data.

info=h5info(fname); % obtain the metadata information

seismic.network=info.Groups(1).Name(2:end); % obtain the network name

lnet=length(seismic.network); % obtain the length of network name

seismic.component=info.Groups(1).Groups(1).Datasets.Name; % obtain the name of the component

seismic.fe=h5read(fname,[info.Groups(2).Name '/fe']); % obtain the sampling frequency of the data

t0_UNIX_timestamp=h5read(fname,[info.Groups(2).Name '/t0_UNIX_timestamp']); % obtain Unix epoch time

seismic.t0=datetime(t0_UNIX_timestamp,'convertfrom','posixtime'); % convert to human readable date and time

ns=length(info.Groups(1).Groups); % obtain the number of stations

nt=info.Groups(1).Groups(1).Datasets.Dataspace.Size; % obtain the number of time samples in a trace

seismic.data=zeros(ns,nt); % initialize the output seismic data

for ii=1:ns
    seismic.sname{ii}=info.Groups(1).Groups(ii).Name(lnet+3:end); % obtain the name of stations
    
    % remove the void location name
    name_temp=seismic.sname{ii};
    if strcmp(name_temp(end-2:end),'.00')
        seismic.sname{ii}=name_temp(1:end-3);
    end
    
    datasetname=['/' seismic.network '/' name_temp '/' seismic.component];
    seismic.data(ii,:)=h5read(fname,datasetname); % obtain seismic data of each station
end


end