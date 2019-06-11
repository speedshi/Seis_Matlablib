function gchkrs(event,trace,tvtsoup,travelp,travels,tcal,show_stk)
% This function is used to display and check the record section for a
% particular input event. The unit of origin time is second.
%
% X-North, Y-East, Z-Vertical down
%
% INPUT:--------------------------------------------------
% event: event origin time and location information, 1*4, Origin_time-X-Y-Z;
% trace: matlab structure, contains seismic data infomation;
% trace.data: seismic data, 2D array, nr*nt;
% trace.dt: scalar, time sampling interval in seconds;
% trace.north: vector, 1*nr, north components of stations;
% trace.east: vector, 1*nr, east compomemts of stations;
% trace.depth: vector, 1*nr, depth components of stations;
% tvtsoup: source imaging positions, 2D array, nn*3, X-Y-Z
% travelp: 2D array, nn*nr, P-wave traveltime table;
% travels: 2D array, nn*nr, S-wave traveltime table;
% tcal: used to calibrate the origin time of the event, in second.

if nargin<6
    tcal=0;
    show_stk=true;
elseif nargin == 6
    show_stk=true;
end

% assemble the station positions, X-Y-Z
recp=[trace.north(:) trace.east(:) trace.depth(:)];

dt=trace.dt; % time sampling interval of the recorded data (s)

% seismic data
skfz=trace.data'; % transpose the data to keep consistent with dispwfscn function


idse=find(abs(tvtsoup(:,1)-event(2))<1e-6 & abs(tvtsoup(:,2)-event(3))<1e-6 & abs(tvtsoup(:,3)-event(4))<1e-6);


soup_mcm=event(2:4); % source position

% record section without t0 calibration
et0=event(1);
dispwfscn(skfz,recp,soup_mcm,dt,et0,travelp(idse,:),travels(idse,:));
title('Record section (without t0 calibration)');

if show_stk
    % show the waveform stack results
    dispwflstk(skfz,dt,et0,travelp(idse,:),3,5);
    title('Stacked waveforms of P-waves (without t0 calibration)');
    
    % show the waveform stack results
    dispwflstk(skfz,dt,et0,travels(idse,:),3,5);
    title('Stacked waveforms of S-waves (without t0 calibration)');
end


% record section with t0 calibration
et0=event(1)+tcal;
dispwfscn(skfz,recp,soup_mcm,dt,et0,travelp(idse,:),travels(idse,:));
title('Record section (with t0 calibration)');

if show_stk
    % show the waveform stack results
    dispwflstk(skfz,dt,et0,travelp(idse,:),3,5);
    title('Stacked waveforms of P-waves (with t0 calibration)');
    
    % show the waveform stack results
    dispwflstk(skfz,dt,et0,travels(idse,:),3,5);
    title('Stacked waveforms of S-waves (with t0 calibration)');
end

end