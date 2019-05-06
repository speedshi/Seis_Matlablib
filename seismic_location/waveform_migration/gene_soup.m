function [soup,snr,ser,sdr,nsnr,nser,nsdr]=gene_soup(north_r,east_r,depth_r,dn,de,dd,precision,fname)
% This function is used to generate source imaging points in 3D space.
%
% Unit: meter.
%
% If null file name ([]) is given for the output file, then the program
% do not output the binary file.
%
% INPUT--------------------------------------------------------------------
% north_r: North component range of source imaing points, vector, 2*1;
% east_r: East component range of source imaging points, vector, 2*1;
% depth_r: Depth component range of source imaging points, vector, 2*1;
% dn: spatial interval of source imaging points in the North direction;
% de: spatial interval of source imaging points in the East direction;
% dd: spatial interval of source imaging points in the Depth direction;
% precision: precision of the output binary file, 'single' or 'double';
% fname: file name of the output binary source position file.
% OUTPUT-------------------------------------------------------------------
% soup: Cartesian coordinates of source imaging points, matrix, nsr*3;
% snr: north coordinates, vector, 1*nsnr;
% ser: east coordinates, vector, 1*nser;
% sdr: depth coordinates, vector, 1*nsdr;
% nsnr: number of imaging points in the north direction, scalar;
% nser: number of imaging points in the east direction, scalar;
% nsdr: number of imaging points in the depth direction, scalar;
% fname: binary source position file.

folder='data'; % name of the folder where output data are stored

% check if the output folder exists, if not, then create it
if ~exist(folder,'dir')
    mkdir(folder);
end

if nargin == 6
    precision='double'; % default output precision
    fname='soupos.dat'; % default output file name
elseif nargin == 7
    fname='soupos.dat'; % default output file name
end

if isempty(precision)
    precision='double';
end

if ~isempty(fname)
    fname=['./' folder '/' fname]; % including the folder
end

% generate source position matrix
snr=north_r(1):dn:north_r(2); % North coordinates
ser=east_r(1):de:east_r(2); % East coordinates
sdr=depth_r(1):dd:depth_r(2); % Depth coordinates

nsnr=max(size(snr)); % number of imaging points in North direction
nser=max(size(ser)); % number of imaging points in East direction
nsdr=max(size(sdr)); % number of imaging points in Depth direction

soup=zeros(nsnr*nser*nsdr,3); % used to store all the source positions, each row is a source position.

% Column 1: North; Column 2: East; Column 3: Depth. 'soup': first change in North axis, then in East axis, then in Depth axis.
for id=1:nsdr
    for ie=1:nser
        for in=1:nsnr
            sid=(id-1)*nsnr*nser+(ie-1)*nsnr+in;
            soup(sid,1)=snr(in); % North coordinates of source imaging points
            soup(sid,2)=ser(ie); % East coordinates of source imaging points
            soup(sid,3)=sdr(id); % Depth coordinates of source imaging points
        end
    end
end

if ~isempty(fname)
    % output binary file
    fid=fopen(fname,'w');
    fwrite(fid,soup,precision);
    fclose(fid);
end

end