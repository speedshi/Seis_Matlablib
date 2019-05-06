function [trace,search,mcm]=mcm_genei(file_seismic,file_stations,file_velocity,search,mcm,precision)
% This function is used to generate the required input files for MCM.
% Unit: meter, m/s, degree.
%
% Accepted seismic data format: h5, SAC
% Accepted station file format: IRIS text
%
% The seismic data of different traces should have the same length, the
% same sampling rate and also start at the same time.
%
% The parameter 'file_seismic' can be a string or cell array which contains
% the file name of the input seismic data.
%
% INPUT--------------------------------------------------------------------
% file_seismic: file name (including path) of the seismic data, a string or cell array;
% file_stations: file name (including path) of the stations, a string;
% file_velocity: file name (including path) of the velocity model, a string;
% search: matlab structure, contains the imaging area information,
% search.north: 1*2, imaging area in the north direction, in meter,
% search.east: 1*2, imaging area in the east direction, in meter,
% search.depth: 1*2, imaging area in the depth direction, in meter;
% search.dn: spatial interval in the north direction, in meter;
% search.de: spatial interval in the east direction, in meter;
% search.dd: spatial interval in the depth direction, in meter;
% mcm: matlab structure, specify MCM parameters, used to generate 'migpara.dat';
% mcm.filter: matlab structure, contains filtering information to filter seismic data;
% mcm.filter.freq: frequency band used to filter the seismic data, a vector containing 1 or 2 elements, in Hz
% mcm.filter.type: filter type, can be 'low', 'bandpass', 'high', 'stop'
% mcm.filter.order: order of Butterworth filter, for bandpass and bandstop designs are of order 2n
% precision: 'single' or 'double', specify the precision of the output
% binary files.
%
% OUTPUT-------------------------------------------------------------------
% waveform.dat: binary file of seismic data for MCM input;
% travelp.dat: binary file of P-wave traveltimes for MCM input;
% travels.dat: binary file of S-wave traveltimes for MCM input;
% soupos.dat: binary file of source imaging positions for MCM input;
% migpara.dat: text file for MCM parameters;
% trace: matlab structure, contain selected data information;
% trace.data: seismic data, 2D array, n_sta*nt;
% trace.dt: time sampling interval of seismic data, in second, scalar;
% trace.name: name of selected stations, vector, 1*n_sta;
% trace.north: north coordinates of selected stations, vector, 1*n_sta;
% trace.east: east coordinates of selected stations, vector, 1*n_sta;
% trace.depth: depth coordinates of selected stations, vector, 1*n_sta;
% trace.t0: matlab datetime, the starting time of traces;
% trace.travelp: P-wave traveltime table, correspond to travelp.dat, 2D array, ns*n_sta;
% trace.travels: S-wave traveltime table, correspond to travels.dat, 2D array, ns*n_sta;
% trace.recp: assambled station positions, 2D array, n_sta*3, N-E-D in meters;
% search.soup: source imaging positions, correspond to soupos.dat, 2D array, ns*3;
% search.snr: north coordinates of imaging points, vector, 1*nsnr;
% search.ser: east coordinates of imaging points, vector, 1*nser;
% search.sdr: depth coordinates of imaging points, vector, 1*nsdr;
% search.nsnr: number of imaging points in the north direction, scalar;
% search.nser: number of imaging points in the east direction, scalar;
% search.nsdr: number of imaging points in the depth direction, scalar;
% mcm: matlab structure, MCM parameters and configuration information.

% set default value
if nargin==4
    mcm=[];
    precision='double';
elseif nargin==5
    precision='double';
end

if isempty(precision)
    precision='double';
end

% set default working directory
if ~isfield(mcm,'workfolder')
    mcm.workfolder='mcm';
end

% check if the working directory exists, if not, then create it
if ~exist(mcm.workfolder,'dir')
    mkdir(mcm.workfolder);
end

% entering the working directory
cd(mcm.workfolder);

folder='./data';
% check if the output folder exists, if not, then create it
if ~exist(folder,'dir')
    mkdir(folder);
end

% process file names of seismic data
if ~isa(file_seismic,'cell')
    % input is characters or string which is the name of a file
    file_seismic=char(file_seismic); % transfer to character vectors
    file_seismic={file_seismic}; % transfer to cell array
end
% read in seismic data
if strcmp(file_seismic{1}(end-2:end),'.h5') || strcmp(file_seismic{1}(end-2:end),'.H5')
    % read in the H5 format data
    seismic=read_seish5(file_seismic);
elseif strcmp(file_seismic{1}(end-3:end),'.sac') || strcmp(file_seismic{1}(end-3:end),'.SAC')
    % read in the SAC format data
    seismic=read_seissac(file_seismic);
else
    error('Unrecognised format of seismic data.');
end

if isfield(mcm,'prefile') && ~isempty(mcm.prefile)
    % use pre-calculated traveltime tables
    load(mcm.prefile,'stations','search'); % load in the 'stations' and the 'search' parameter
    
    % obtain MCM required input files
    % generate binary file of source imaging positions
    [search.soup,search.snr,search.ser,search.sdr,search.nsnr,search.nser,search.nsdr]=gene_soup(search.north,search.east,search.depth,search.dn,search.de,search.dd,precision);
    
else
    % need to calculate traveltime tables
    % read in station information
    stations=read_stations(file_stations); % read in station information in IRIS text format
        
    % read in velocity infomation
    model=read_velocity(file_velocity); % read in velocity model, now only accept homogeneous and layered model
        
    % obtain MCM required input files
    % generate binary file of source imaging positions
    [search.soup,search.snr,search.ser,search.sdr,search.nsnr,search.nser,search.nsdr]=gene_soup(search.north,search.east,search.depth,search.dn,search.de,search.dd,precision);
    
    % generate traveltime tables
    stations=gene_traveltime(model,stations,search,precision,[],[]);
end


if ~isfield(mcm,'filter')
    mcm.filter=[];
end

% generate binary files for seismic data and traveltimes
trace=gene_wavetime(seismic,stations,mcm.filter,precision);

% assemble the station positions, X-Y-Z i.e. North-East-Depth
trace.recp=[trace.north(:) trace.east(:) trace.depth(:)];

% generate text files for mcm parameters
if isfield(mcm,'migtp')
    mcm.nre=size(trace.data,1); % number of stations
    mcm.nsr=size(search.soup,1); % number of imaging points
    mcm.dfname='waveform.dat'; % file name of seismic data
    mcm.dt=trace.dt; % time sampling interval, in second
    mcm.tdatal=(size(trace.data,2)-1)*trace.dt; % time length of the whole seismic data in second (s)
    mcm.vthrd=0.001; % threshold value for identifying seismic event in the migration volume
    mcm.spaclim=0; % the space limit in searching for potential seismic events, in meter (m)
    mcm.timelim=0; % the time limit in searching for potential seismic events, in second (s)
    mcm.nssot=1; % the maximum number of potential seismic events can be accept for a single origin time
    gene_migpara(mcm); % generate the text file
end


addrs=sprintf('%s/info.mat',folder); % name of output file
% output matlab format data, can be used for later
save(addrs,'trace','search','mcm');

if ~isfield(mcm,'run')
    mcm.run=0;
end

% check if need to run the MCM program
switch mcm.run
    case 1
        fprintf('Run MCM Fortran-OpenMP program.\n');
        
    case 2
        fprintf('Run MCM Matlab program.\n');
        
    case 3
        fprintf('Run MCM Matlab test version program.\n');
        
        if ~isfield(mcm.test,'timerg') || isempty(mcm.test.timerg)
            % set default time range
            mcm.test.timerg=[mcm.test.t0; mcm.test.t0+seconds((size(trace.data,2)-1)*trace.dt)];
        end
        
        migv=runmcm_matlab_test(trace,mcm,search,mcm.test);
    otherwise
        fprintf('No MCM program is running.\n');
end

cd('..');

end