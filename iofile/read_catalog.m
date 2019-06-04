function catalog=read_catalog(file_cata,timerg,search)
% This function is used to read in the catalog information. The catalog
% file is in IRIS text format.
%
% Use SI unit: meter, degree. Elevation or depth is relative to the
% sea-level. For depth, + means below the sea-level, - means above the
% sea-level.
%
% INPUT--------------------------------------------------------------------
% file_cata: catalog file name;
% timerg: matlab datetime, vetor: 2*1, used to select events in the time
% range;
% search: structure, used to select events in the specific zone;
% search.north: 1*2, imaging area in the north direction, in meter;
% search.east: 1*2, imaging area in the east direction, in meter;
% search.depth: 1*2, imaging area in the depth direction, in meter;
%
% OUTPUT-------------------------------------------------------------------
% catalog: matlab structure which contains corresponding catalog infomation
% catalog.time: origin time of earthquakes, in matlab datetime format;
% catalog.latitude: latitude in degree;
% catalog.longitude: longitude in degree;
% catalog.elevation: elevation in meter;
% catalog.north: North Cartisian coordinate in meter;
% catalog.east: East Cartisian coordinate in meter;
% catalog.depth: Depth Cartisian coordinate in meter;
% catalog.magnitude: magnitude of earthquakes;

% set default parameters
if nargin == 1
    timerg=[];
    search=[];
elseif nargin == 2
    search=[];
end

opts=detectImportOptions(file_cata);
opts.VariableNames{1}='EventID'; % set the name of the first colume
catainfo=readtable(file_cata,opts); % read in all the information from file


catalog.latitude=catainfo.Latitude; % obtain latitude info
catalog.longitude=catainfo.Longitude; % obtain longitude info
catalog.elevation=catainfo.Depth_Km*-1000; % obtain elevation info, note the unit transfer

% obtain Cartisian coordinates using default wgs84Ellipsoid system
[catalog.east,catalog.north,catalog.depth]=geod2cart(catalog.latitude,catalog.longitude,catalog.elevation);

% obtain magnitude info
catalog.magnitude=catainfo.Magnitude;

% obtain origin time info
catalog.time=datetime(catainfo.Time,'Inputformat','yyyy-MM-dd''T''HH:mm:ss.SSSSSS');

% select events within correct time range
if ~isempty(timerg)
    indx = catalog.time>=timerg(1) & catalog.time<=timerg(2);
    catalog.time=catalog.time(indx);
    catalog.latitude=catalog.latitude(indx);
    catalog.longitude=catalog.longitude(indx);
    catalog.elevation=catalog.elevation(indx);
    catalog.north=catalog.north(indx);
    catalog.east=catalog.east(indx);
    catalog.depth=catalog.depth(indx);
    catalog.magnitude=catalog.magnitude(indx);
end

% select events within specific zone
if ~isempty(search)
    indx = catalog.north>=search.north(1) & catalog.north<=search.north(2) & catalog.east>=search.east(1) & catalog.east<=search.east(2) & catalog.depth>=search.depth(1) & catalog.depth<=search.depth(2);
    catalog.time=catalog.time(indx);
    catalog.latitude=catalog.latitude(indx);
    catalog.longitude=catalog.longitude(indx);
    catalog.elevation=catalog.elevation(indx);
    catalog.north=catalog.north(indx);
    catalog.east=catalog.east(indx);
    catalog.depth=catalog.depth(indx);
    catalog.magnitude=catalog.magnitude(indx);
end

end