function output_catatxt1(catalog,fname,para)
% This function is used to output the catalog in text format.
% The output file can then be used for GMT plotting.
% NOTE: the output unit is exactly the same as the input unit, no transfer here.
%
% INPUT--------------------------------------------------------------------
%       catalog: structure, contains catalog information, e.g.
%       catalog.time, catalog.magnitude, catalog.latitude,
%       catalog.longitude, catalog.depth, catalog.distance...
%       fname: string, filename for the output file;
%       para: structure, contains parameters to format the output file;
%


% set default parameters
if nargin < 2
    fname='catalog.txt';
    para.rftime='toordinal';
    para.distance=false;
end

if ~isfield(para,'rftime')
    para.rftime='toordinal';
end

if ~isfield(para,'distance')
    para.distance=false;
end

n=length(catalog.time); % number of earthquakes in the input catalog



if strcmp(para.rftime,'toordinal')
    % relative days count from the date AD 01/01/01, i.e.
    % equivalent to python datetime.toordinal() function.
    rtimes=datenum(catalog.time)-366;
else
    % use days relative to a input datetime
    rtimes=days(catalog.time-para.rftime);
end



fid=fopen(fname,'wt');

% output the text file
if para.distance
    % output distances, the catalog must contain distance information
    for ii=1:n
        % note the format of the output file:
        % longitude-latitude-depth-origin times-magnitude-distances
        fprintf(fid,'%f %f %f %f %f %f \n',catalog.longitude(ii),catalog.latitude(ii),catalog.depth(ii),rtimes(ii),catalog.magnitude(ii),catalog.distance(ii));
    end
else
    % do not output distances
    for ii=1:n
        % note the format of the output file: longitude-latitude-depth-origin times-magnitude
        fprintf(fid,'%f %f %f %f %f \n',catalog.longitude(ii),catalog.latitude(ii),catalog.depth(ii),rtimes(ii),catalog.magnitude(ii));
    end
end

fclose(fid);

end