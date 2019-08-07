function ttime=waveforms_show(event,trace,search,tcal)
% This function is used to show the corresponding waveforms of an
% earthquake.
%
% Input distance is in meter.
%
% INPUT--------------------------------------------------------------------
% event: vector, 1*4, Origin_time-North-East-Depth, origin time is relative
% to the trace.t0 in second;
% trace: matlab structure, contain selected data information;
% trace.data: seismic data, 2D array, n_sta*nt;
% trace.dt: time sampling interval of seismic data, in second, scalar;
% trace.name: name of selected stations, vector, 1*n_sta;
% trace.north: north coordinates of selected stations, vector, 1*n_sta;
% trace.east: east coordinates of selected stations, vector, 1*n_sta;
% trace.depth: depth coordinates of selected stations, vector, 1*n_sta;
% trace.t0: matlab datetime, the starting time of traces;
% trace.travelp: P-wave traveltime table, 2D array, ns*n_sta;
% trace.travels: S-wave traveltime table, 2D array, ns*n_sta;
% trace.recp: assambled station positions, 2D array, n_sta*3, N-E-D in meters;
% search: matlab structure, contains the imaging area information;
% search.soup: source imaging positions, correspond to soupos.dat, 2D array, ns*3;
% tcal: scalar to calibrate the origin time of the event, in second;
%
% OUTPUT-------------------------------------------------------------------
% ttime: structure, time information of seismic signals at different
% stations (traveltimes and arrivaltimes) expressed relative to the
% starting time of seismic data; 
% Note arrival-times are the calibrated arrival-times;
% ttime.travelp: P-wave travel-times in second, shape: n_sta*1;
% ttime.travels: S-wave travel-times in second, shape: n_sta*1;
% ttime.arritp: P-wave arrival-times in second, shape: n_sta*1;
% ttime.arrits: S-wave arrival-times in second, shape: n_sta*1;


% set default parameters
if nargin<4
    tcal=0;
end


% fine the index of source in the traveltime table
idxe=abs(search.soup(:,1)-event(2))<1e-6 & abs(search.soup(:,2)-event(3))<1e-6 & abs(search.soup(:,3)-event(4))<1e-6;

% obtain the P-wave traveltimes for this event location
travelp=trace.travelp(idxe,:);

% obtain the S-wave traveltimes for this event location
travels=trace.travels(idxe,:);

% obtain seismic data
data=trace.data';

% assemble the station positions, North-East-Depth, transfer from m to km
recp=[trace.north(:) trace.east(:) trace.depth(:)]/1000;

% assemble the source positions, North-East-Depth, transfer from m to km
soup=event(2:4)/1000;

% obtain time sampling interval
dt=trace.dt;

% obtain origin time of the seismic event
stos=event(1);

% obtain the name of stations
name=trace.name;

% obtain the starting time of seismic data
t0=trace.t0;

% obtain record section without t0 calibration
dispwfscn_mn(data,recp,soup,dt,stos,travelp,travels,name,t0);
title('Record section (without t0 calibration)');

% obtain origin time of the seismic event after t0 calibration
stos=event(1)+tcal;

% obtain record section with t0 calibration
dispwfscn_mn(data,recp,soup,dt,stos,travelp,travels,name,t0);
title('Record section (with t0 calibration)');


% sort out the output travel-time information
ttime.travelp=travelp;
ttime.travels=travels;

ttime.arritp=event(1)+tcal+travelp;
ttime.arrits=event(1)+tcal+travels;


end