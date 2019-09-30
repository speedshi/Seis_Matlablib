function [sname,nsta]=read_staname(file_seismic)
% This file is used to read the trace/station names in the seismic data
% file
%
% INPUT--------------------------------------------------------------------
% file_seismic: string, name of the input seismic file;
%
% OUTPUT-------------------------------------------------------------------
% sname: cell array, the names of the stations/traces in the seismic data;
% nsta: scalar, the total number of stations/traces;



% obtain the station names in the HDF5 file--------------------------------
info=h5info(file_seismic); % read the information of the HDF5 file

num_net=length(info.Groups(:))-1; % the number of networks, note the last one in the Groups should be '_metadata', thus is not a network

% obtain the total number of stations in the HDF5 file
nsta=0;
for inet=1:num_net
    nsta=nsta+length(info.Groups(inet).Groups);
end

% initial cell array
network=cell(nsta,1); % network name
sname=cell(nsta,1); % station/trace name

id_tra=0; % trace/station index
for inet=1:num_net
    % Read in seismic data for each network
    ns=length(info.Groups(inet).Groups); % obtain the number of stations in the current network
    
    for ii=1:ns
        
        id_tra=id_tra+1; % the current trace index
        
        network{id_tra}=info.Groups(inet).Name(2:end); % obtain the current network name
        lnet=length(network{id_tra}); % obtain the length of the network name
        sname{id_tra}=info.Groups(inet).Groups(ii).Name(lnet+3:end); % obtain the name of stations
        
        % remove the void location name
        name_temp=sname{id_tra};
        if strcmp(name_temp(end-2:end),'.00')
            sname{id_tra}=name_temp(1:end-3);
        end
        
    end
    
end


end