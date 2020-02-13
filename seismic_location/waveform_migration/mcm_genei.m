function [trace,search,mcm]=mcm_genei(file,search,mcm,precision)
% This function is used to generate the required input files for MCM.
% Unit: meter, m/s, degree.
%
% Accepted seismic data format: HDF5, SAC, MAT;
% Accepted station file format: IRIS text;
% Accepted catalog file format: IRIS text;
%
% The seismic data of different traces should have the same length, the
% same sampling rate and also start at the same time.
%
% The parameter 'file.seismic' can be a string or cell array which contains
% the file name of the input seismic data.
%
% INPUT--------------------------------------------------------------------
% file: matlab structure, contains the file names of the input data;
% file.seismic: file name (including path) of the seismic data, a string or
%               cell array; Or if file.seismic is a structure array, then it
%               contains all needed information about seismic data, no need
%               to load from file.
% file.stations: file name (including path) of the stations, a string; Or
%                if file.stations is a structure array, then it contains
%                all needed information about seismic stations, no need to
%                load from file.
% file.velocity: file name (including path) of the velocity model, a
%                string; Or if file.velocity is a structure array, then it
%                contains all needed information about velocity model, no
%                need to load from file.
% search: matlab structure, contains the imaging area information;
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
% mcm.utmstruct: struture, the UTM parameter for coordinate transfermation;
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


% set default parameters---------------------------------------------------
if nargin==2
    mcm=[];
    precision='double';
elseif nargin==3
    precision='double';
end

if isempty(precision)
    precision='double';
end

% set default seismic data component to perform migration
if ~isfield(mcm,'component')
    mcm.component='Z';
end

% set default working directory
if ~isfield(mcm,'workfolder')
    mcm.workfolder='mcm';
end

% set utmstruct
if ~isfield(mcm,'utmstruct')
    mcm.utmstruct = [];
end


% check if the working directory exists, if not, then create it
mcm.workfolder=['./' mcm.workfolder];
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

% read the seismic data
if isstruct(file.seismic)
    % input is a structure array, contain everything about seismic data, no
    % need to load from file
    seismic = file.seismic;
else
    % input is the file name of seismic data, need to load information from
    % file
    seismic=read_seis(file.seismic);
end

% check if need to reset the t0 of the seismic data
if isfield(mcm,'datat0')
    seismic.t0=mcm.datat0;
end

% show the spectrogram of the original seismic data
if isfield(mcm,'stationid') && ~isempty(mcm.stationid)
    ispectrogram(seismic,mcm.stationid);
end

if isfield(mcm,'prefile') && ~isempty(mcm.prefile)
    % use pre-calculated traveltime tables and source imaging points
    load(mcm.prefile,'stations','search'); % load in the 'stations' and the 'search' parameter
    
else
    % need to calculate traveltime tables
    % read in station information
    if isstruct(file.stations)
        % input is a structure array, contain everything about station
        % information, no need to load from file
        stations = file.stations;
    else
        % input is the file name of station data, need to load information
        % from file
        stations=read_stations(file.stations,[],mcm.utmstruct); % read in station information in IRIS text format
        mcm.utmstruct = stations.utmstruct; % reset the utmstruct, since this is the first point to transform coordinates
    end
    
    % read in velocity infomation
    if isstruct(file.velocity)
        % input is a structure array, contain everything about velocity
        % model, no need to load from file
        model = file.velocity;
    else
        % input is the file name of velocity model, need to load
        % information from file
        model=read_velocity(file.velocity); % read in velocity model, now only accept homogeneous and layered model
    end
    
    % obtain MCM required input files
    % generate binary file of source imaging positions
    [search.soup,search.snr,search.ser,search.sdr,search.nsnr,search.nser,search.nsdr]=gene_soup(search.north,search.east,search.depth,search.dn,search.de,search.dd,precision,[]);
    
    % generate traveltime tables
    stations=gene_traveltime(model,stations,search,precision,[],[]);
end


if ~isfield(mcm,'filter')
    mcm.filter=[];
end

% generate binary files for seismic data and traveltimes
trace=gene_wavetime(seismic,stations,mcm.filter,precision,[],[],[]);

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
end


if ~isfield(mcm,'run')
    mcm.run=-1;
end


mcm.datat0=seismic.t0; % the starting time (t0) of seismic data, datetime

% mcm.dtimerg: time range of the seismic data
mcm.dtimerg=[mcm.datat0; mcm.datat0+seconds((size(trace.data,2)-1)*trace.dt)];


% check if need to run the MCM program
switch mcm.run
    case -1
        fprintf('No MCM program is running. Just generate the input files for MCM Fortran program.\n');
        gene_mcmifiles(trace,search,mcm,precision);
        
    case 0
        fprintf('Run MCM parameter testing program at the source location.\n');
        
        % obtain the information of the specified earthquake
        earthquake=get_earthquake(mcm,search);
        
        % obtain the migration trace at the input earthquake location
        [s_pro,n_var]=mcm_test_para(trace,mcm,search,earthquake);
        
    case 1
        fprintf('Run MCM frequency band testing.\n');
        
        % obtain the information of the specified earthquake
        earthquake=get_earthquake(mcm,search);
        
        % search and show results for different frequency band
        [s_pro,n_var]=mcm_test_freqband(file.seismic,stations,mcm,earthquake,search);
        
    case 2
        fprintf('Run MCM testing program with the catalog input.\n');
        
        % obtain the information of the specified earthquake
        earthquake=get_earthquake(mcm,search);
        
        % obtain the searching origin time serials
        mcm.st0=(earthquake.t0-mcm.test.twind):mcm.dt0:(earthquake.t0+mcm.test.twind);
        
        % run mcm around the time range of the input earthquake
        migv=runmcm_matlab_test(trace,mcm,search,earthquake);
        
    case 3
        fprintf('Run MCM Matlab testing program with an input time range.\n');
        
        if isfield(mcm,'test') && isfield(mcm.test,'mtrg')
            % chech if the input time range are feasible
            if mcm.test.mtrg(1)<mcm.dtimerg(1) || mcm.test.mtrg(1)>mcm.test.mtrg(2) || mcm.test.mtrg(2)>=mcm.dtimerg(2)
                error('Incorrect input for mcm.test.mtrg! Not applicable.');
            end
            
            % obtain the searching origin time serials
            time_start=seconds(mcm.test.mtrg(1)-mcm.datat0); % starting time for migration, relative to data_t0 in second
            time_end=seconds(mcm.test.mtrg(2)-mcm.datat0); % ending time for migration, relative to data_t0 in second
            mcm.st0=time_start:mcm.dt0:time_end;
        else
            % set default searching time range
            mcm=detmst0(mcm,trace);
        end
        
        % run mcm in the input time range
        migv=runmcm_matlab_test(trace,mcm,search);
        
    case 4
        fprintf('Run MCM Matlab program.\n');
        
        [migv,mcm]=runmcm_matlab(trace,mcm,search);
        
    case 99
        fprintf('Run MCM Fortran-OpenMP program.\n');
        
    otherwise
        fprintf('No MCM program is running. Just load in data.\n');
end

cd('..');

end