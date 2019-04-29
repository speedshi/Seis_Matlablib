function stations=read_stations(fname)
% This function is used to read the station file to obtain station
% information such as the name and location.
%
% The depth parameter in the station file represent:
% The local depth or overburden of the instrument's location. For downhole
% instruments, the depth of the instrument under the surface ground level.
% For underground vaults, the distance from the instrument to the local
% ground level above.
%
% Location is expressed in Cartesian coordinate system. The convention of
% coordinate systerm follows the definition in Aki & Richards 2002, Fig
% 4.20 and BOX 4.4. Right-hand coordinate systerm in Cartesian coordinate,
% with X-North, Y-East, Z-Vertical down (Depth). Unit: meter.
% For the axis of Vertical down (depth), the sea-level is 0, and above the
% sea-level is negative, and below the sea-level is positive; can be
% interpret as depth relative to the sea-level.
%
% The station file uses IRIS text format. The first line is the header
% which explains the meaning of each column thereafter.
% INPUT-----------------------------------------------
% fname: file name including path of the station file;
% OUTPUT----------------------------------------------
% stations: structure, contains  metadata of each station;
% stations.name: cell array, the name of each station;
% stations.north: vector, North components of the position of each station;
% stations.east: vector, East components of the position of each station;
% stations.depth: vector, Depth components of the position of each station.

opts=detectImportOptions(fname);
opts.VariableNames{1}='Network'; % set the name of the first colume
stall=readtable(fname,opts); % read in all the information
nnr=size(stall,1); % number of all station items in the file

% name of the station
snamelist{1}=stall.Station{1};
Location=stall.Location{1};
if isempty(Location)
    Location='00';
end

% stations.name{1}=[snamelist{1} '.' Location];
stations.name{1}=snamelist{1};

nr=1;

% Obtain Cartesian coordinates of the station
[stations.north(1),stations.east(1),stations.depth(1)]=geod2cart(stall.Latitude(1),stall.Longitude(1),stall.Elevation(1)-stall.Depth(1));

for ii=2:nnr
    
    if ~ismember(stall.Station{ii},snamelist)
        
        nr=nr+1; % total number of different stations
        
        snamelist{nr}=stall.Station{ii};
        Location=stall.Location{ii};
        if isempty(Location)
            Location='00';
        end
        %stations.name{nr}=[snamelist{nr} '.' Location];
        stations.name{nr}=snamelist{nr};
        
        % Obtain Cartesian coordinates of the stations
        [stations.north(nr),stations.east(nr),stations.depth(nr)]=geod2cart(stall.Latitude(ii),stall.Longitude(ii),stall.Elevation(ii)-stall.Depth(ii));
    end
    
end



end