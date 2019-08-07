function stations=gene_traveltime(model,stations,search,precision,fname_p,fname_s)
% This function is used to generate traveltime tables and binary files which
% are required by MCM.
%
% Unit for input parameters are: meters, m/s.
%
% The depth parameter is relative to the sea-level. Positive: belew the
% sea-level; Negative: above the sea-level.
%
% If null file names ([]) are given for the output files, then the program
% do not output the binary files.
%
% INPUT--------------------------------------------------------------------
% model: matlab structure, contains velocity and layering information;
% model.rsd0: reference starting depth of the model, i.e. depth of free surface;
% model.thickness: vector, thickness of each layer;
% model.vp: vector, P-wave velocities of each layer;
% model.vs: vector, S-wave velocities of each layer.
% stations: matlab structure, contains station position information;
% stations.name: cell array, 1*nre, the name of each station;
% stations.north: vector, 1*nre, North components of the position of each station;
% stations.east: vector, 1*nre, East components of the position of each station;
% stations.depth: vector, 1*nre, Depth components of the position of each station.
% search: matlab structure, contains the imaging area information;
% search.soup: Cartesian coordinates of source imaging points (X-Y-Z or N-E-D), matrix, nsr*3;
% precision: strings, 'single' or 'double', specifiy the outout presicion;
% fname_p: output binary file name for P-wave traveltimes;
% fname_s: output binary file name for S-wave traveltimes.
%
% OUTPUT-------------------------------------------------------------------
% stations: matlab structure;
% stations.travelp: P-wave traveltime table, in second, 2D array, nsr*nre;
% stations.travels: S-wave traveltime table, in second, 2D array, nsr*nre;
% stations_traveltime_search.mat: matlat binary file, contains 'stations' 'search';
% fname_p: binary file of P-wave traveltimes for MCM input;
% fname_s: binary file of S-wave traveltimes for MCM input.

folder='./data'; % name of the folder where output data are stored

% check if the output folder exists, if not, then create it
if ~exist(folder,'dir')
    mkdir(folder);
end

% set default names of the output binary files
if nargin==3
    precision='double';
    fname_p='travelp.dat';
    fname_s='travels.dat';
elseif nargin==4
    fname_p='travelp.dat';
    fname_s='travels.dat';
end

if isempty(precision)
    precision='double';
end

if ~isempty(fname_p)
    fname_p=[folder '/' fname_p]; % including the folder
end

if ~isempty(fname_s)
    fname_s=[folder '/' fname_s]; % including the folder
end

stations.travelp=[];
stations.travels=[];

soup=search.soup; % positions of source imaging points

% correct for the non-zero reference starting depth of the model,
% because by default the 'tvtcalrt' program assumes the free surface of the
% model is at 0 depth (sea-level).
recp_z=stations.depth-model.rsd0; % correct depth values for stations
soup(:,3)=search.soup(:,3)-model.rsd0; % correct depth values for source imaging points

recp=[stations.north(:) stations.east(:) recp_z(:)]; % N-E-D

if length(model.thickness)==1
    % homogeneous model
    [stations.travelp,stations.travels,~,~]=tvtcalrt_homo(model.vp/1000.,model.vs/1000.,soup/1000.,recp/1000.);
else
    % layered model
    if ~isempty(model.vp)
        [stations.travelp,~,~]=tvtcalrt_ly(model.vp/1000.,model.thickness/1000.,soup/1000.,recp/1000.); % note the unit transfer
    end
    if ~isempty(model.vs)
        [stations.travels,~,~]=tvtcalrt_ly(model.vs/1000.,model.thickness/1000.,soup/1000.,recp/1000.); % note the unit transfer
    end
end

% output mat format files
save('stations_traveltime_search.mat','stations','search');

% output binary files
if ~isempty(stations.travelp) && ~isempty(fname_p)
    % for P-wave traveltimes
    fid=fopen(fname_p,'w');
    fwrite(fid,stations.travelp,precision);
    fclose(fid);
end

if ~isempty(stations.travels) && ~isempty(fname_s)
    % for S-wave traveltimes
    fid=fopen(fname_s,'w');
    fwrite(fid,stations.travels,precision);
    fclose(fid);
end

end