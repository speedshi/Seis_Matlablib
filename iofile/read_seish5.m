function seismic=read_seish5(fname,component)
% This function is used to read seismic data of HDF5 format.
%
% The seismic data of different traces and components must have the same
% length, the same sampling rate and also start at the same time.
%
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
% component: cell array, can be 'Z' 'N' 'E', component of seismic data to read;
% OUTPUT-------------------------------------------------
% seismic: structure, contains seismic data and metadata of each station;
% seismic.network: cell array, the name of the networks;
% seismic.component: cell array, the name of the data component (usually N, E or Z);
% seismic.name: cell array, 1*ns, contains the name of each station;
% seismic.fe: scaler, the sampling frequency of the data, in Hz;
% seismic.dt: scaler, the time sample interval of the data, in second;
% seismic.t0: matlab datetime, the origin time of the seismic data;
% seismic.zdata: 2D array, ns*nt, seismic data of Z component;
% seismic.ndata: 2D array, ns*nt, seismic data of N component;
% seismic.edata: 2D array, ns*nt, seismic data of E component;



% set default parameters
if nargin<2
    component={'Z'}; % read in Z component by default
end

fname=char(fname); % file name of HDF5 file

info=h5info(fname); % read the information of the HDF5 file

seismic.component=component; % set which component to read


% obtain metadata information----------------------------------------------
metaname='/_metadata'; % i.e. "info.Groups(end).Name"
seismic.fe=h5read(fname,[metaname '/fe']); % obtain the sampling frequency of the data
seismic.dt=1.0/seismic.fe; % obtain time sample interval of the data, in second
t0_UNIX_timestamp=h5read(fname,[metaname '/t0_UNIX_timestamp']); % obtain Unix epoch time
seismic.t0=datetime(t0_UNIX_timestamp,'convertfrom','posixtime'); % convert to human readable date and time


% Read in seismic data-----------------------------------------------------
num_net=length(info.Groups(:))-1; % the number of networks, note the last one in the Groups should be '_metadata', thus is not a network

% obtain the total number of stations in the HDF5 file
nsta=0;
for inet=1:num_net
    nsta=nsta+length(info.Groups(inet).Groups);
end

id_compo=nan;
for ii=1:length(info.Groups(1).Groups(1).Datasets)
    if info.Groups(1).Groups(1).Datasets(ii).Name == seismic.component{1}
        id_compo=ii;
    end
end

if isnan(id_compo)
    error('Component: %s not include in the data file.',seismic.component{1});
end

nt=info.Groups(1).Groups(1).Datasets(id_compo).Dataspace.Size; % obtain the number of time samples in a trace

num_comp=length(seismic.component); % obtain the number of components to read in
for ic=1:num_comp
    % initialize the output seismic data
    if strcmp(seismic.component{ic},'Z') || strcmp(seismic.component{ic},'z')
        % for Z component
        seismic.zdata=zeros(nsta,nt);
    elseif strcmp(seismic.component{ic},'N') || strcmp(seismic.component{ic},'n')
        % for N component
        seismic.ndata=zeros(nsta,nt);
    elseif strcmp(seismic.component{ic},'E') || strcmp(seismic.component{ic},'e')
        % for E component
        seismic.edata=zeros(nsta,nt);
    else
        error('Unrecognised input component: %s for reading seismic data!', seismic.component{ic});
    end
end

id_tra=0; % trace index
for inet=1:num_net
    % Read in seismic data for each network
    ns=length(info.Groups(inet).Groups); % obtain the number of stations in the current network
    
    for ii=1:ns
        
        id_tra=id_tra+1; % the current trace index
        
        seismic.network{id_tra}=info.Groups(inet).Name(2:end); % obtain the current network name
        lnet=length(seismic.network{id_tra}); % obtain the length of the network name
        seismic.name{id_tra}=info.Groups(inet).Groups(ii).Name(lnet+3:end); % obtain the name of stations
        
        % remove the void location name
        name_temp=seismic.name{id_tra};
        if strcmp(name_temp(end-2:end),'.00')
            seismic.name{id_tra}=name_temp(1:end-3);
        end
        
        % read in seismic data for input components
        for ic=1:num_comp
            % obtain the path to the seismic data
            datasetname=['/' seismic.network{id_tra} '/' name_temp '/' seismic.component{ic}];
            
            % read in data
            if strcmp(seismic.component{ic},'Z') || strcmp(seismic.component{ic},'z')
                % for Z component
                seismic.zdata(id_tra,:)=h5read(fname,datasetname); % obtain seismic data of each station
            elseif strcmp(seismic.component{ic},'N') || strcmp(seismic.component{ic},'n')
                % for N component
                seismic.ndata(id_tra,:)=h5read(fname,datasetname); % obtain seismic data of each station
            elseif strcmp(seismic.component{ic},'E') || strcmp(seismic.component{ic},'e')
                % for E component
                seismic.edata(id_tra,:)=h5read(fname,datasetname); % obtain seismic data of each station
            else
                error('Unrecognised input component: %s for reading seismic data!', seismic.component{ic});
            end
            
        end
        
    end
    
end

end