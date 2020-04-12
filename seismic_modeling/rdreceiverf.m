function station=rdreceiverf(fname)
% Read the receiver file.
% Coordinate convention: X-North, Y-East, Z-Depth (vertical down);
% Unit: meter.
%
% INPUT:-------------------------------------------------
% fname: file name of the receiver file.
% OUTPUT:------------------------------------------------
% station: structure, contains information about stations;
% station.nr: scalar, number of stations;
% station.recp: matrix: nr*3, X-Y-Z positions of stations;
% station.name: cell array: nr*1, name of stations;
% stations.north: vector, North components of the position of each station;
% stations.east: vector, East components of the position of each station;
% stations.depth: vector, Depth components of the position of each station;


if nargin<1
    fname='receiver.dat';
end

sdata=importdata(fname,' ',1);
station.nr=size(sdata.data,1); % number of receivers
station.recp=sdata.data; % X-Y-Z locations
station.name=sdata.textdata(2:end); % station names
station.north =sdata.data(:,1); % North coordinates
station.east=sdata.data(:,2); % East coordinates
station.depth=sdata.data(:,3); % Depth coordinates

end