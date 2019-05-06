function sltsname=sltordotpsta(file_stations,point,num,file_out)
% This function is used to select, order and output the name of stations.
% Stations is ordered with respect to the distance to a certain point.
%
% INPUT--------------------------------------------------------------------
% file_stations: the name of station file, in IRIS text format;
% point: matlab structure, the position of reference point;
% point.latitude, in degree; point.longitude, in degree;
% or point.north, in meter; point.east, in meter;
% num: the number of stations to select, order and output;
% file_out: the output file name.
% OUTPUT:------------------------------------------------------------------
% sltsname: slected station name list in correct order;
% file_out: text file which contains the name of the stations.

% set default parameters
if nargin == 3
    file_out='./station_select.txt';
elseif nargin == 2
    file_out='./station_select.txt';
    num=0;
end

% read in station information
stations=read_stations(file_stations); % read in station information in IRIS text format

% transfer to Cartisian coordinates
if isfield(point,'latitude') && ~isfield(point,'north') 
    [point.east,point.north,~]=geod2cart(point.latitude,point.longitude,0);
end

% calculate the horizontal distance between the stations and the input point
stations.hdistance=sqrt((stations.north-point.north).^2+(stations.east-point.east).^2);

% sort the horizontal distance in ascending order
[~,indx]=sort(stations.hdistance);

% select stations which are closer to the point
if num>0
    % select the first 'num' stations according to the distance
    sltsname=stations.name(indx(1:num));
else
    % save all the stations in the input station file
    sltsname=stations.name(indx);
end

nn=length(sltsname);
% output the selected station names
fid=fopen(file_out,'wt');
for ii=1:nn
    fprintf(fid,sltsname{ii});
    fprintf(fid,'\n');
end
fclose(fid);

end