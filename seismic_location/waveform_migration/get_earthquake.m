function earthquake=get_earthquake(mcm,search)
% This function is used to obtain the specified earthquake information from
% the catalog.
%
% INPUT--------------------------------------------------------------------
% mcm: structure, containing the catalog and earthquake information;
% mcm.test.cataname: string, the file name of the input catalog, in IRIS text format;
% mcm.dtimerg: datetime, 2*1, the time range of the seismic data;
% mcm.test.cataid: the event id of the earthquake in the catalog;
% mcm.datat0: datetime, the starting time (t0) of seismic data;
% search: structure, used to select events in the specific zone;
% search.north: 1*2, imaging area in the north direction, in meter;
% search.east: 1*2, imaging area in the east direction, in meter;
% search.depth: 1*2, imaging area in the depth direction, in meter;
%
% OUTPUT-------------------------------------------------------------------
% earthquake: structure, containing the earthquake information;
% earthquake.north: north component of earthquake location, in meter;
% earthquake.east: east component of earthquake location, in meter;
% earthquake.depth: depth component of earthquake location, in meter;
% earthquake.t0: relative origin time of the earthquake, in second, relative
% to the starting time of the seismic data;

if nargin == 1
    search=[];
end


% read in catalog data
catalog=read_catalog(mcm.test.cataname,mcm.dtimerg,search);

% obtain the relative tested origin time (relative to data t0), in second
earthquake.t0=seconds(catalog.time(mcm.test.cataid)-mcm.datat0);

% obtain the location of the tested earthquake, in meter
earthquake.north=catalog.north(mcm.test.cataid); % north component
earthquake.east=catalog.east(mcm.test.cataid); % east component
earthquake.depth=catalog.depth(mcm.test.cataid); % depth component


end