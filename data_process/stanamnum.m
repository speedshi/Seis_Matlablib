function [sname,snum]=stanamnum(file_seismic,file_stations)
% This function is used to obtain the number and names of stations in the data.
%
% INPUT--------------------------------------------------------------------
% file_seismic: file name of the H5 seismic file;
% file_stations: file name of the station file in IRIS text format;
%
% OUTPUT-------------------------------------------------------------------
% sname: name of the stations;
% snum: number of the stations;

% set default value
if nargin==1
    file_stations=[];
end

info=h5info(file_seismic); % obtain the metadata information of seismic data

seismic.network=info.Groups(1).Name(2:end); % obtain the network name

lnet=length(seismic.network); % obtain the length of network name

ns=length(info.Groups(1).Groups); % obtain the number of stations

for ii=1:ns
    seismic.sname{ii}=info.Groups(1).Groups(ii).Name(lnet+3:end); % obtain the name of stations
    
    % remove the void location name
    name_temp=seismic.sname{ii};
    if strcmp(name_temp(end-2:end),'.00')
        seismic.sname{ii}=name_temp(1:end-3);
    end
end

if isempty(file_stations)
    % no input station file
    sname=seismic.sname;
    snum=ns;
else
    % have input station file
    stations=read_stations(file_stations);
    snum=0;
    for ii=1:ns
        if ismember(seismic.sname{ii},stations.name)
            snum=snum+1;
            sname{snum}=seismic.sname{ii};
        end
    end
    
end

% print information
fprintf('Find %d stations.\n',snum);
fprintf('Station names:');
for ii=1:snum
    fprintf(' %s',sname{ii});
end
fprintf('.\n');

end