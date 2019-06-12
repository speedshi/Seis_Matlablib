function seismic=read_seis(fname)
% This function is used to call different functions to read in the seismic
% data of different formats.
%
% The function utilize the suffix to decide the format of the input seismic
% data and then call the corresponding functions to read the seismic data.
%
% Accept seismic data format: HDF5, SAC.
%
% INPUT--------------------------------------------------------------------
% fname: file name of the input seismic data, can be a string or cell
% array which contains the names of all seismic files;
%
% OUTPUT-------------------------------------------------------------------
% seismic: structure, contains seismic data and metadata of each station;
% seismic.name: cell array, 1*ns, contains the name of each station;
% seismic.fe: scaler, the sampling frequency of the data, in Hz;
% seismic.dt: scaler, the time sample interval of the data, in second;
% seismic.t0: matlab datetime, the origin time of the seismic data;
% seismic.data: 2D array, ns*nt, contains seismic data;
% seismic.network: string, the name of the network;
% seismic.component: character, the name of the data component (usually N, E or Z);


% process file names of seismic data
if ~isa(fname,'cell')
    % input is characters or string which is the name of a file
    fname=char(fname); % transfer to character vectors
    fname={fname}; % transfer to cell array
end


% read in seismic data
if strcmp(fname{1}(end-2:end),'.h5') || strcmp(fname{1}(end-2:end),'.H5')
    % read in the H5 format data
    seismic=read_seish5(fname);
elseif strcmp(fname{1}(end-3:end),'.sac') || strcmp(fname{1}(end-3:end),'.SAC')
    % read in the SAC format data
    seismic=read_seissac(fname);
else
    error('Unrecognised format of seismic data.');
end



end