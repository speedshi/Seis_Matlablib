function seismic=read_seismat(fname)
% This function is used to read seismic data of MAT format.
%
%
% For the MAT format data, it should contain:
% seismic data, must contain at least one of the following;
% zdata: 2D array, ns*nt, continuous seismic data, (must/optional);
% ndata: 2D array, ns*nt, continuous seismic data, (must/optional);
% edata: 2D array, ns*nt, continuous seismic data, (must/optional);
% dt: scaler, the time sample interval of the data, in second, (must);
% name: cell array, 1*ns, contains the name of each station;
% fe: scaler, the sampling frequency of the data, in Hz;
% t0: matlab datetime, the origin time of the seismic data;
% network: cell array, the name of the networks;
%
%
% INPUT--------------------------------------------------------------------
% fname: path and file name of the input seismic data;
% para: structure, controlling parameters;
% para.component: specify to load which component of seismic data;
% OUTPUT-------------------------------------------------------------------
% seismic: structure, contains seismic data and metadata of each station;
% seismic.network: cell array, 1*ns, the name of the networks;
% seismic.name: cell array, 1*ns, contains the name of each station;
% seismic.fe: scaler, the sampling frequency of the data, in Hz;
% seismic.dt: scaler, the time sample interval of the data, in second;
% seismic.t0: matlab datetime, the origin time of the seismic data;
% seismic.zdata: 2D array, ns*nt, continuous seismic data of Z component;
% seismic.ndata: 2D array, ns*nt, continuous seismic data of N component;
% seismic.edata: 2D array, ns*nt, continuous seismic data of E component;
%--------------------------------------------------------------------------


% load data
seismic=load(fname);


% chech field and set default one if there is no input---------------------
if ~isfield(seismic,'zdata') && ~isfield(seismic,'ndata') && ~isfield(seismic,'edata')
    error('The MAT format is incorrect! At least one of the zdata/ndata/edata need to exis.');
elseif isfield(seismic,'zdata')
    ns=size(seismic.zdata,1); % number of stations
elseif isfield(seismic,'ndata')
    ns=size(seismic.ndata,1); % number of stations
elseif isfield(seismic,'edata')
    ns=size(seismic.edata,1); % number of stations
end


if ~isfield(seismic,'dt')
    error('The MAT format is incorrect! No: dt.')
end

if ~isfield(seismic,'fe')
    seismic.fe=1.0/seismic.dt;
end


if ~isfield(seismic,'name')
    for ir=1:ns
        seismic.name{ir}=['S' num2str(ir)]; % default name is just 'S + the trace number'
    end
end

if ~isfield(seismic,'network')
    seismic.network=cell(1,ns); % default network is empty cell array
end

if ~isfield(seismic,'t0')
    seismic.t0=0; % default starting time
end



end