function mcm=detmst0(mcm,trace)
% This function is used to determine the searching origin times for the
% migrating.
% Time in second.
%
% INPUT--------------------------------------------------------------------
% mcm: structure, contains the MCM parameters;
% mcm.tdatal: time length of the whole seismic data in second (s);
% mcm.dt0: time sampling interval of searching origin times in second;
% mcm.tpwind: P-phase time window length in second, scalar;
% mcm.tswind: S-phase time window length in second, scalar;
% mcm.phasetp: specify seismic phase used for migration, scalar;
% trace: matlab structure, contains seismic data information;
% trace.data: seismic data, 2D array, nre*nt;
% trace.dt: time sampling interval, in second;
% trace.travelp: P-wave traveltime table, 2D array, nsr*nre;
% trace.travels: S-wave traveltime table, 2D array, nsr*nre;
%
% OUTPUT-------------------------------------------------------------------
% mcm.st0: searched origin times of MCM, in second (relative to start time
% of input seismic data), vector, nst0*1;


% determine the searched origin times of MCM; in second (relative to start
% time of input seismic data)
time_start=0;



if mcm.phasetp==0
    % only use P-phase
    time_end=mcm.tdatal-max(trace.travelp(:))-mcm.tpwind;
elseif mcm.phasetp==1
    % only use S-phase
    time_end=mcm.tdatal-max(trace.travels(:))-mcm.tswind;
elseif mcm.phasetp==2
    % use both P- and S-phase
    maxwind=max(mcm.tpwind,mcm.tswind);
    time_end=mcm.tdatal-max(trace.travels(:))-maxwind;
else
    error('Incorrect input for mcm.phasetp!');
end


mcm.st0=time_start:mcm.dt0:time_end;


end