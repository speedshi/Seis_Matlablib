function wrtasgeo(soup,recpx,recpy,fname)
% This function is used to write the geometry of surface array and source
% into file 'asgeometry'. The 'asgeopetry' is used as input file for
% generating full waveform data later.
% Usage: 
% input:
%         soup---position of source points(X-Y-Z) (m), could be more than
%         one source, each row is a source;
%         recpx--X position range of surface arrays (m);
%         recpy--Y position range of surface arrays (m), must have the same
%         length as 'recpx'.
% output: (saved at '../data/')
%           a text file whose name is 'asgeometry'.
% The file format of the 'asgeometry' file is:
%  row 1: X position of surface geophones
%  row 2: Y position of surface geophones
%  row 3: X position of source points
%  row 4: Y position of source points
%  row 5: Z position of source points
% The unit of position in file 'asgeometry' is km.

if nargin<4
    fname='asgeometry.txt';
end

% trandfer unit: m -> km
tsoup=soup/1000;
trax=recpx/1000;
tray=recpy/1000;

fid1=fopen(fname,'wt');

% extend position of surface array
[recpxg,recpyg]=meshgrid(trax,tray);
recx=recpxg(:); recy=recpyg(:);

fprintf(fid1,'%.4f \t',recx); % write X position of surface geophones
fprintf(fid1,'\n'); % new line
fprintf(fid1,'%.4f \t',recy); % write Y position of surface geophones
fprintf(fid1,'\n'); % new line
fprintf(fid1,'%.4f \t',tsoup(:,1)); % write X position of source points
fprintf(fid1,'\n'); % new line
fprintf(fid1,'%.4f \t',tsoup(:,2)); % write Y position of source points
fprintf(fid1,'\n'); % new line
fprintf(fid1,'%.4f \t',tsoup(:,3)); % write Z position of source points

fclose(fid1);

end