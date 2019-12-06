function output_catatxt1(catalog,fname,para)
% This function is used to output the catalog in text format.
% The output file can then be used for GMT plotting.
%
% INPUT--------------------------------------------------------------------
%       catalog: structure, contains catalog information, e.g.
%       catalog.time, catalog.magnitude, catalog.latitude,
%       catalog.longitude, catalog.depth...
%       fname: string, filename for the output file;
%       para: structure, contains parameters to format the output file;
%


% set default parameters
if nargin < 2
    fname='catalog.txt';
    para=[];
end


n=length(catalog.time); % number of earthquakes in the input catalog

if ~isempty(para)
    if isfield(para,'rftime')
        if strcmp(para.rftime,'toordinal')
            % relative days count from the date AD 01/01/01, i.e.
            % equivalent to python datetime.toordinal() function.
            rtimes=datenum(catalog.time)-366;
        else
            % use days relative to a input datetime
            rtimes=days(catalog.time-para.rftime);
        end
    end
    
end


fid=fopen(fname,'wt');

% output the text file
for ii=1:n
    % note the format of the output file: longitude-latitude-depth-origin
    % times-magnitude
    fprintf(fid,'%f %f %f %f %f \n',catalog.longitude(ii),catalog.latitude(ii),catalog.depth(ii),rtimes(ii),catalog.magnitude(ii));
end

fclose(fid);

end