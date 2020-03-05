function stations=read_stations(fname,select,utmstruct)
% This function is used to read the station file to obtain station
% information such as the name and location.
%
% The Depth parameter in the IRIS station file represent:
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
%
% If there are Cartesian coordinates in the file (North, East, Depth), then
% do not read and transfer the geographic coordinates, read in and use the
% Cartesian coordinates directly.
%
% INPUT-----------------------------------------------
% fname: file name including path of the station file;
% select: structure, select stations that fullfill the requirements;
% select.name: only select the named stations;
% utmstruct: the UTM parameters for coordinate transfermation;
%
% OUTPUT----------------------------------------------
% stations: structure, contains  metadata of each station;
% stations.name: cell array, the name of each station;
% stations.north: vector, North components of the position of each station;
% stations.east: vector, East components of the position of each station;
% stations.depth: vector, Depth components of the position of each station;
% stations.latitude: vector, latitude in degree of each station;
% stations.longitude: vector, longitude in degree of each station;
% stations.elevation: vector, elevation in meter of each station;
% stations.utmstruct: the UTM parameters for coordinate transfermation;


if nargin < 2
    select=[];
    utmstruct=[];
elseif nargin < 3
    utmstruct=[];
end

opts=detectImportOptions(fname);
opts.VariableNames{1}='Network'; % set the name of the first colume
stall=readtable(fname,opts); % read in all the information
nnr=size(stall,1); % number of all station items in the file

% name of the station
snamelist{1}=stall.Station{1};
% Location=stall.Location{1};
% if isempty(Location)
%     Location='00';
% end

% stations.name{1}=[snamelist{1} '.' Location];
stations.name{1}=snamelist{1};

if ismember('Latitude',stall.Properties.VariableNames)
    % Obtain geographic coordinates of the station
    stations.latitude(1)=stall.Latitude(1);
    stations.longitude(1)=stall.Longitude(1);
    stations.elevation(1)=stall.Elevation(1)-stall.Depth(1);
else
    % no geographic coordinates, then the Cartesian coordinates must exist
    stations.north(1)=stall.North(1);
    stations.east(1)=stall.East(1);
    stations.depth(1)=stall.Depth(1);
end


nr=1;
for ii=2:nnr
    
    if ~ismember(stall.Station{ii},snamelist)
        
        nr=nr+1; % total number of different stations
        
        snamelist{nr}=stall.Station{ii};
%         Location=stall.Location{ii};
%         if isempty(Location)
%             Location='00';
%         end
        %stations.name{nr}=[snamelist{nr} '.' Location];
        stations.name{nr}=snamelist{nr};
        
        if ismember('Latitude',stall.Properties.VariableNames)
            % Obtain the geographic information of the stations
            stations.latitude(nr)=stall.Latitude(ii);
            stations.longitude(nr)=stall.Longitude(ii);
            stations.elevation(nr)=stall.Elevation(ii)-stall.Depth(ii);
        else
            % no geographic coordinates, then the Cartesian coordinates must exist
            stations.north(nr)=stall.North(ii);
            stations.east(nr)=stall.East(ii);
            stations.depth(nr)=stall.Depth(ii);
        end
    end
    
end


% Obtain Cartesian coordinates of the stations
% Note, in order to use the same UTM zone to convert coordinate system, all
% stations must be transfered together!
if ismember('Latitude',stall.Properties.VariableNames)
    [stations.east,stations.north,stations.depth,stations.utmstruct] = ...
    geod2cart(stations.latitude,stations.longitude,stations.elevation,utmstruct);
end


if ~isempty(select)
    % select the stations that fullfill the requirements
    lindx=ismember(stations.name,select.name);
    
    stations.name=stations.name(lindx);
    stations.east=stations.east(lindx);
    stations.north=stations.north(lindx);
    stations.depth=stations.depth(lindx);
    
    if ismember('Latitude',stall.Properties.VariableNames)
        stations.latitude=stations.latitude(lindx);
        stations.longitude=stations.longitude(lindx);
        stations.elevation=stations.elevation(lindx);
    end
    
end


end