function migv=runmcm_matlab_test(trace,mcm,search,earthquake)
% This function is used to set up and run MCM Matlab test version.
% X-North  Y-East  Z-Depth(Vertical down)
%
% INPUT--------------------------------------------------------------------
% trace: matlab structure, contains seismic data information;
% trace.data: seismic data, 2D array, nre*nt;
% trace.dt: time sampling interval, in second;
% trace.travelp: P-wave traveltime table, 2D array, nsr*nre;
% trace.travels: S-wave traveltime table, 2D array, nsr*nre;
% trace.recp: assambled station positions, 2D array, nre*3, N-E-D in meters;
% trace.name: name of selected stations, vector, 1*n_sta;
% trace.t0: matlab datetime, the starting time of traces;
% mcm: matlab structure, contains MCM parameters;
% mcm.phasetp: specify seismic phase used for migration, scalar;
% mcm.tpwind: P-phase time window length in second, scalar;
% mcm.tswind: S-phase time window length in second, scalar;
% mcm.st0: searched origin times of MCM, in second (relative to start time
% of input seismic data), vector, nst0*1;
% mcm.dt0: time sampling interval of searching origin times in second;
% search: matlab structure, describe the imaging area;
% search.north: 1*2, imaging area in the north direction, in meter,
% search.east: 1*2, imaging area in the east direction, in meter,
% search.depth: 1*2, imaging area in the depth direction, in meter;
% search.dn: spatial interval in the north direction, in meter;
% search.de: spatial interval in the east direction, in meter;
% search.dd: spatial interval in the depth direction, in meter;
% search.soup: source imaging positions, 2D array, nsr*3, in meter;
% search.nsnr: number of imaging points in the north direction, scalar;
% search.nser: number of imaging points in the east direction, scalar;
% search.nsdr: number of imaging points in the depth direction, scalar;
% earthquake: matlab structure, contains the location and origin time of the earthquake;
% earthquake.north: scalar, earthquake location in north direction, in meter;
% earthquake.east: scalar, earthquake location in east direction, in meter;
% earthquake.depth: scalar, earthquake location in depth direction, in meter;
% earthquake.t0: scalar, relative earthquake origin time, in second,
% relative to the origin time of the seismic data;
%
% OUTPUT-------------------------------------------------------------------
% migv: migration volume, 4D array, shape: nsnr*nser*nsdr*nst0.


% set default parameters
if nargin<4
    earthquake=[];
end

% calculate the characteristic function
trace.data=transpose(cal_charfunc(trace.data',mcm.cfuntp));

% run the MCM test
if mcm.phasetp==2
    % use both P and S phase
    fprintf('Use P-phase of %f s window and S-phase of %f s window for MCM.\n',mcm.tpwind,mcm.tswind);
    migv=wave_migration_kernel(trace,mcm,search);
elseif mcm.phasetp==0
    % use only P-phase
    fprintf('Use P-phase of %f s window for MCM.\n',mcm.tpwind);
    mcm.txwind=mcm.tpwind;
    trace.travelx=trace.travelp;
    migv=wave_migration_kernel_x(trace,mcm,search);
elseif mcm.phasetp==1
    % use only S-phase
    fprintf('Use S-phase of %f s window for MCM.\n',mcm.tswind);
    mcm.txwind=mcm.tswind;
    trace.travelx=trace.travels;
    migv=wave_migration_kernel_x(trace,mcm,search);
else
    error('Incorrect input for mcm.phasetp, only accept: 0, 1, 2.');
end

% show the migration result
show_mcmmigres(migv,search,trace,mcm,earthquake);

end