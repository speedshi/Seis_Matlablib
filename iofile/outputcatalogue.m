function outputcatalogue(events,rtime,ut_zone,fname,stname,tvtsoup,travelp,travels)
% This function is used to output event location result to a formated
% catalogue file.
% INPUT:-------------------------------------------------
% events: event location information, nsr*4, time-X-Y-Z. Time is relative to
% the starting time of recorded data in second, XYZ coordinates is in km
% and in a specific UTM zones;
% rtime: the referenced starting time, struct.
% ut_zone: the UTM zone of this area, string;
% fname: the name of the output file;
% stname: the name of the stations, must keep consistent with travelp and travels;
% tvtsoup: potential source positions, nps*3, X-Y-Z, in km;
% travelp: travel-time table of P-wave, nps*nre, in second;
% travels: travel-time table of S-wave, nps*nre, in second.

nsr=size(events,1); % number of seismic events
nre=length(stname); % number of stations

% output README file which describes the format of each output file.
fid0=fopen('README','wt');
fprintf(fid0,'"_cata": seismic catalogue file.\r\n');
fprintf(fid0,'    Column 1: Seismic event id;\r\n');
fprintf(fid0,'    Column 2: Origin time UTC;\r\n');
fprintf(fid0,'    Column 3: Relative origin time, expressed in seconds after %4d-%02d-%02dT%02d:%02d:%07.4f UTC;\r\n',rtime.year,rtime.month,rtime.day,rtime.hour,rtime.minute,rtime.second);
fprintf(fid0,'    Column 4: Latitude (in degrees);\r\n');
fprintf(fid0,'    Column 5: Longitude (in degrees);\r\n');
fprintf(fid0,'    Column 6: Depth (in km), relative to the sea level, positive value means below the sea level;\r\n');
fprintf(fid0,'    Column 7: UTM North coordinate (in km);\r\n');
fprintf(fid0,'    Column 8: UTM East coordinate (in km);\r\n');
fprintf(fid0,'\r\n');
fprintf(fid0,'"_ptime": Calculated P-wave arrival times at different stations.\r\n');
fprintf(fid0,'    Station list: '); fprintf(fid0,'%s ',stname{:}); fprintf(fid0,'\r\n');
fprintf(fid0,'    Column 1: Seismic event id;\r\n');
fprintf(fid0,'    Column 2: Relative arrival time at station 1 (%s), expressed in seconds after %4d-%02d-%02dT%02d:%02d:%07.4f UTC;\r\n',stname{1},rtime.year,rtime.month,rtime.day,rtime.hour,rtime.minute,rtime.second);
fprintf(fid0,'    Column 3: Relative arrival time at station 2 (%s), expressed in seconds after %4d-%02d-%02dT%02d:%02d:%07.4f UTC;\r\n',stname{2},rtime.year,rtime.month,rtime.day,rtime.hour,rtime.minute,rtime.second);
fprintf(fid0,'    ...\r\n');
fprintf(fid0,'\r\n');
fprintf(fid0,'"_stime": Calculated S-wave arrival times at different stations.\r\n');
fprintf(fid0,'    Station list: '); fprintf(fid0,'%s ',stname{:}); fprintf(fid0,'\r\n');
fprintf(fid0,'    Column 1: Seismic event id;\r\n');
fprintf(fid0,'    Column 2: Relative arrival time at station 1 (%s), expressed in seconds after %4d-%02d-%02dT%02d:%02d:%07.4f UTC;\r\n',stname{1},rtime.year,rtime.month,rtime.day,rtime.hour,rtime.minute,rtime.second);
fprintf(fid0,'    Column 3: Relative arrival time at station 2 (%s), expressed in seconds after %4d-%02d-%02dT%02d:%02d:%07.4f UTC;\r\n',stname{2},rtime.year,rtime.month,rtime.day,rtime.hour,rtime.minute,rtime.second);
fprintf(fid0,'    ...\r\n');
fclose(fid0);
%------------------------------------------------------------


% output event catalogue-------------------------------
fid1=fopen([fname '_cata'],'wt');
fprintf(fid1,'Event_id  UTC origin time           Relative origin time(s)    Latitude     Longitude      Depth(km)     North(km)        East(km)\r\n');
% define the map projection structure
utmstruct = defaultm('utm');
utmstruct.zone = ut_zone;
utmstruct.geoid = wgs84Ellipsoid;
utmstruct = defaultm(utmstruct);
% set the reference UTC time
rt_utc=datetime(rtime.year,rtime.month,rtime.day,rtime.hour,rtime.minute,rtime.second);
for iis=1:nsr
    % find the correct UTC origin times
    ot_utc=datevec(datetime(events(iis,1),'ConvertFrom','epochtime','Epoch',rt_utc)); % orgin time vector
    % transfer from map to geographic coordinates
    [latitude,longitude]=minvtran(utmstruct,events(iis,3)*1000,events(iis,2)*1000); % note need to transfer from km to m for the input X-Y coordinates.
    fprintf(fid1,'%5d     %4d-%02d-%02dT%02d:%02d:%07.4f     %11.4f          %10.4f    %10.4f    %10.4f    %12.4f    %12.4f\r\n',iis,ot_utc(1),ot_utc(2),ot_utc(3),ot_utc(4),ot_utc(5),ot_utc(6),events(iis,1),latitude,longitude,events(iis,4),events(iis,2),events(iis,3));
end
fclose(fid1);
%------------------------------------------------------------


% output P-wave arrival times-------------------------
fid2=fopen([fname '_ptime'],'wt');
fprintf(fid2,'event_id      ');
for iir=1:nre
    fprintf(fid2,'%s         ',stname{iir});
end
fprintf(fid2,'\r\n');
for iis=1:nsr
    idse=abs(tvtsoup(:,1)-events(iis,2))<1e-6 & abs(tvtsoup(:,2)-events(iis,3))<1e-6 & abs(tvtsoup(:,3)-events(iis,4))<1e-6; % logical indexing
    art=events(iis,1)+travelp(idse,:); % arrival time = origin time + traveltime
    fprintf(fid2,'%5d     ',iis);
    for iir=1:nre
        fprintf(fid2,'%11.4f  ',art(iir));
    end
    fprintf(fid2,'\r\n');
end
fclose(fid2);
%------------------------------------------------------------


% output S-wave arrival times-------------------------
fid3=fopen([fname '_stime'],'wt');
fprintf(fid3,'event_id      ');
for iir=1:nre
    fprintf(fid3,'%s         ',stname{iir});
end
fprintf(fid3,'\r\n');
for iis=1:nsr
    idse=abs(tvtsoup(:,1)-events(iis,2))<1e-6 & abs(tvtsoup(:,2)-events(iis,3))<1e-6 & abs(tvtsoup(:,3)-events(iis,4))<1e-6; % logical indexing
    art=events(iis,1)+travels(idse,:); % arrival time = origin time + traveltime
    fprintf(fid3,'%5d     ',iis);
    for iir=1:nre
        fprintf(fid3,'%11.4f  ',art(iir));
    end
    fprintf(fid3,'\r\n');
end
fclose(fid3);
%------------------------------------------------------------
end