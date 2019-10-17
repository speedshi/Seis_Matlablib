function seismic=read_seismat(fname)
% This function is used to read seismic data of MAT format.
%
%
% For the MAT format data, it should contain:
% data: 2D array, ns*nt, continuous seismic data, (must);
% dt: scaler, the time sample interval of the data, in second, (must);
% name: cell array, 1*ns, contains the name of each station;
% fe: scaler, the sampling frequency of the data, in Hz;
% t0: matlab datetime, the origin time of the seismic data;
% network: cell array, the name of the networks;
% component: cell array, the name of the data component (usually N, E or Z);
%
%
% INPUT--------------------------------------------------
% fname: path and file name of the input seismic data;
% OUTPUT-------------------------------------------------
% seismic: structure, contains seismic data and metadata of each station;
% seismic.network: cell array, 1*ns, the name of the networks;
% seismic.component: cell array, the name of the data component (usually N, E or Z);
% seismic.name: cell array, 1*ns, contains the name of each station;
% seismic.fe: scaler, the sampling frequency of the data, in Hz;
% seismic.dt: scaler, the time sample interval of the data, in second;
% seismic.t0: matlab datetime, the origin time of the seismic data;
% seismic.data: 2D array, ns*nt, continuous seismic data;


% load data
seismic=load(fname);

% chech field and set default one if there is no input---------------------
if ~isfield(seismic,'data')
    error('The MAT format is incorrect! No: data.')
end

if ~isfield(seismic,'dt')
    error('The MAT format is incorrect! No: dt.')
end

if ~isfield(seismic,'fe')
    seismic.fe=1.0/seismic.dt;
end

ns=size(seismic.data,1); % number of stations
if ~isfield(seismic,'name')
   for ir=1:ns
      seismic.name{ir}=['S' num2str(ir)]; % default name is just 'S + the trace number'
   end
end

if ~isfield(seismic,'network')
    seismic.network=cell(1,ns); % default network is empty cell array
end

if ~isfield(seismic,'t0')
    seismic.t0=datetime('2015-01-01 00:00:00'); % default starting time
end

if ~isfield(seismic,'component')
    seismic.component={'Z'}; % default component
end


end