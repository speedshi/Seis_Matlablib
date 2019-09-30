function [sname,snum]=stanamnum(file_seismic,file_stations)
% This function is used to obtain the number and names of stations in the HDF5 file.
% The stations must exist in the 'file_stations' as well.
%
% INPUT--------------------------------------------------------------------
% file_seismic: file name of the H5 seismic file;
% file_stations: file name of the station file in IRIS text format;
%
% OUTPUT-------------------------------------------------------------------
% sname: name of the stations;
% snum: total number of the stations;


% set default value
if nargin==1
    file_stations=[];
end



% read in the names of stations/trace in the seismic file
[seismic.name, nsta]=read_staname(file_seismic);



% check if all the stations are in the input 'file_station' file-----------
if isempty(file_stations)
    % no input station file
    sname=seismic.name;
    snum=nsta;
else
    % have input station file
    stations=read_stations(file_stations);
    snum=0;
    for ii=1:nsta
        if ismember(seismic.name{ii},stations.name)
            snum=snum+1;
            sname{snum}=seismic.name{ii};
        end
    end
    
end

% print information
fprintf('Find %d stations.\n',snum);
if snum>0
    % find stations
    fprintf('Station names:');
    for ii=1:snum
        fprintf(' %s',sname{ii});
    end
    fprintf('.\n');
else
    % does not find any stations
    sname=[];
end

end