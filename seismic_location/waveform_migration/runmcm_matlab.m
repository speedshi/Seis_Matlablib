function [migv,mcm]=runmcm_matlab(trace,mcm,search)
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
% mcm: matlab structure, contains MCM parameters;
% mcm.phasetp: specify seismic phase used for migration, scalar;
% mcm.tpwind: P-phase time window length in second, scalar;
% mcm.tswind: S-phase time window length in second, scalar;
% mcm.dt0: time sampling interval of searching origin times in second;
% search: matlab structure, describe the imaging area;
% search.soup: source imaging positions, 2D array, nsr*3, in meter;
%
% OUTPUT-------------------------------------------------------------------
% migv: migration volume, 2D array, shape: nst0*5; for each row:
% Origin_time-North-East-Depth-Migration_value;
% mcm: strucure, contains the MCM parameters;


% determine the searched origin times
mcm=detmst0(mcm,trace);

% calculate the characteristic function
trace.data=transpose(cal_charfunc(trace.data',mcm.cfuntp));

% run the MCM program
if mcm.phasetp==2
    % use both P and S phase
    fprintf('Use P-phase of %f s window and S-phase of %f s window for MCM.\n',mcm.tpwind,mcm.tswind);
    migv=wave_migration_kernel_r(trace,mcm,search);
elseif mcm.phasetp==0
    % use only P-phase
    fprintf('Use P-phase of %f s window for MCM.\n',mcm.tpwind);
    mcm.txwind=mcm.tpwind;
    trace.travelx=trace.travelp;
    migv=wave_migration_kernel_x_r(trace,mcm,search);
elseif mcm.phasetp==1
    % use only S-phase
    fprintf('Use S-phase of %f s window for MCM.\n',mcm.tswind);
    mcm.txwind=mcm.tswind;
    trace.travelx=trace.travels;
    migv=wave_migration_kernel_x_r(trace,mcm,search);
else
    error('Incorrect input for mcm.phasetp, only accept: 0, 1, 2.');
end




end