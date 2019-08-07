function gene_mcmifiles(trace,search,mcm,precision,fname_r,fname_d,fname_p,fname_s)
% This function is used to generate the correct input files for MCM Fortran
% program.
%
% If the file names are set to be empty, the the binary file will not be
% output.
%
% INPUT--------------------------------------------------------------------
% trace: matlab structure, contain selected data information;
% trace.data: seismic data, 2D array, n_sta*nt;
% trace.travelp: P-wave traveltime table, 2D array, ns*n_sta;
% trace.travels: S-wave traveltime table, 2D array, ns*n_sta;
% search: matlab structure, contains the imaging area information;
% search.soup: source imaging positions, correspond to soupos.dat, 2D array, ns*3;
% mcm: matlab structure, specify MCM parameters, used to generate 'migpara.dat';
% precision: 'single' or 'double', the precision of the output files;
% fname_r: characters, the file name of the binary source position file;
% fname_d: characters, the file name of the binary seismic data file;
% fname_p: characters, the file name of the binary P-wave traveltime file;
% fname_s: characters, the file name of the binary S-wave traveltime file;
%
% OUTPUT-------------------------------------------------------------------
% soupos.dat: binary file of source imaging positions;
% waveform.dat: binary file of seismic data;
% travelp.dat: binary file of P-wave traveltime table;
% travels.dat: binary file of S-wave traveltime table;
% migpara.dat: text file of MCM parameters;
% info.mat: matlab format data of seismic data and the related information;



% set default parameters
if nargin == 3
    precision='double';
    fname_r='soupos.dat';
    fname_d='waveform.dat';
    fname_p='travelp.dat';
    fname_s='travels.dat';
elseif nargin == 4
    fname_r='soupos.dat';
    fname_d='waveform.dat';
    fname_p='travelp.dat';
    fname_s='travels.dat';
elseif nargin == 5
    fname_d='waveform.dat';
    fname_p='travelp.dat';
    fname_s='travels.dat';
elseif nargin == 6
    fname_p='travelp.dat';
    fname_s='travels.dat';
elseif nargin == 7
    fname_s='travels.dat';
end

if isempty(precision)
    precision='double';
end


folder='./data'; % name of the folder where output data are stored

% check if the output folder exists, if not, then create it
if ~exist(folder,'dir')
    mkdir(folder);
end


% obtain MCM required input files
% generate binary file of source imaging positions
if ~isempty(fname_r)
    fid=fopen([folder '/' fname_r],'w');
    fwrite(fid,search.soup,precision);
    fclose(fid);
end


% generate binary file of seismic data
if ~isempty(fname_d) && ~isempty(trace.data)
    fid=fopen([folder '/' fname_d],'w');
    fwrite(fid,trace.data,precision);
    fclose(fid);
end


% generate binary file of P-wave traveltime table
if isfield(trace,'travelp') && ~isempty(trace.travelp)
    fid=fopen([folder '/' fname_p],'w');
    fwrite(fid,trace.travelp,precision);
    fclose(fid);
end


% generate binary file of S-wave traveltime table
if isfield(trace,'travels') && ~isempty(trace.travels)
    fid=fopen([folder '/' fname_s],'w');
    fwrite(fid,trace.travels,precision);
    fclose(fid);
end


% output text file of MCM parameters
gene_migpara(mcm); % generate the text file


% output matlab format data of seismic data and the related information, can be used for later
addrs=sprintf('%s/info.mat',folder); % name of output file
save(addrs,'trace','search','mcm');


end