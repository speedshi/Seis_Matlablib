function event = show_mcmmigres(migv,search,trace,mcm,earthquake)
% This function is used to show the migration results of MCM.
%
% INPUT--------------------------------------------------------------------
% migv: migration volume, 4D array, shape: nsnr*nser*nsdr*nst0;
% search: matlab structure, containing the imaging area information;
% search.north: 1*2, imaging area in the north direction, in meter,
% search.east: 1*2, imaging area in the east direction, in meter,
% search.depth: 1*2, imaging area in the depth direction, in meter;
% search.soup: source imaging positions, 2D array, ns*3, in meter;
% search.nsnr: number of imaging points in the north direction, scalar;
% search.nser: number of imaging points in the east direction, scalar;
% search.nsdr: number of imaging points in the depth direction, scalar;
% trace: matlab structure, contains seismic data information;
% trace.data: seismic data, 2D array, nre*nt;
% trace.dt: time sampling interval, in second;
% trace.travelp: P-wave traveltime table, 2D array, nsr*nre;
% trace.travels: S-wave traveltime table, 2D array, nsr*nre;
% trace.recp: assambled station positions, 2D array, n_sta*3, N-E-D in meters;
% trace.name: name of selected stations, vector, 1*nre;
% trace.t0: matlab datetime, the starting time of traces;
% mcm: structure, contains parameters for mcm;
% mcm.st0: vector, nst0*1, searched origin times of MCM, in second
% (relative to the start time of input seismic data);
% mcm.phasetp: specify seismic phase used for migration, scalar;
% mcm.tpwind: time window in second for P-phase;
% mcm.tdatal: time length of the whole seismic data in second (s);
% mcm.utmstruct: struture, the UTM parameter for coordinate transfermation;
% earthquake: matlab structure, contains the location and origin time of
% the earthquake, this input can have null input;
% earthquake.north: scalar, earthquake location in north direction, in meter;
% earthquake.east: scalar, earthquake location in east direction, in meter;
% earthquake.depth: scalar, earthquake location in depth direction, in meter;
% earthquake.t0: scalar, relative earthquake origin time, in second,
% relative to the origin time of the seismic data;
%
% OUTPUT-------------------------------------------------------------------
% event: structure, contains information about the located event;
% event.t0: origin time of the event, in datetime format;
% event.latitude: latitude of the event, in degree;
% event.longitude: longitude of the event, in degree;
% event.depth: depth of the event, in meter;


% set default parameters
if nargin<5
    earthquake=[];
end


% show and compare migration result

wfmstk_c=permute(migv,[4 1 2 3]); % note here exchange dimensions

if ~isempty(earthquake)
    % have earthquake input
    soup_cata=[earthquake.north earthquake.east earthquake.depth]; % location of the earthquake, in meter
else
    % no earthquake input
    soup_cata=[];
end

taxis = trace.t0 + seconds(mcm.st0); % searched origin times (absolute times)
[tn,xn,yn,zn]=migmaxplt(wfmstk_c,soup_cata/1000,search.north/1000,search.east/1000,search.depth/1000,taxis); % note the unit transfer, m->km

idse=sub2ind([search.nsnr search.nser search.nsdr],xn,yn,zn); % location index for the MCM

event.t0 = trace.t0 + seconds(mcm.st0(tn));  % origin time of the located seismic event (datatime format)
time_str = datestr(event.t0,'yyyy-mm-dd HH:MM:SS.FFF');

[event.latitude,event.longitude]=minvtran(mcm.utmstruct,search.soup(idse,2),search.soup(idse,1)); % transfer Cartesian (m) to Geodetic coordinate (degree)
event.depth = search.soup(idse,3);

fprintf('Maximum coherence value: %f.\n',max(wfmstk_c(:))); % maximum coherency value in the volume
fprintf('Origin time: %s.\n',time_str); % origin time of the located event
fprintf('Event location: X-%f, Y-%f, Z-%f m.\n',search.soup(idse,:)); % located event locations
fprintf('Latitude: %f; Longitude: %f; Depth: %f m.\n',event.latitude,event.longitude,event.depth); % located event locations


lwin=30; % left window length for showing seismic data, in second
rwin=60; % right window length for showing seismic data, in second

nre=size(trace.data,1); % number of stations used in MCM



% Display results, for MCM without origin tima calibration
% display the arrival times of seismic event on the recorded seismic data
et0=mcm.st0(tn); % the located origin time

% make sure the set time range is not out of boundary
if lwin>et0
    lwinc=et0;
else
    lwinc=lwin;
end
if rwin>mcm.tdatal-et0
    rwinc=mcm.tdatal-et0;
else
    rwinc=rwin;
end

net0r=round((et0-lwinc)/trace.dt+1):round((et0+rwinc)/trace.dt+1); % sample index for extracting waveforms from continuous data
exwfm_t0 = trace.t0 + seconds(et0-lwinc); % t0 time of the extracted data, in datetime format (absolute times)
exwfm=transpose(trace.data(:,net0r)); % extracted waveforms
extmx = trace.t0 + seconds((net0r-1)*trace.dt);  % time axis for waveforms (absolute times)
for ire=1:nre
    figure; plot(extmx,exwfm(:,ire),'k'); hold on;
    if mcm.phasetp==0 || mcm.phasetp==2
        % show P-phase time window
        ttmin=trace.t0+seconds(trace.travelp(idse,ire)+et0);
        ttmax=ttmin+seconds(mcm.tpwind);
        vmin=min(exwfm(:,ire));
        vmax=max(exwfm(:,ire));
        fill([ttmin ttmax ttmax ttmin],[vmin vmin vmax vmax],'b','linestyle','none','FaceAlpha',0.35); hold on;
    end
    if mcm.phasetp==1 || mcm.phasetp==2
        % show S-phase time window
        ttmin=trace.t0+seconds(trace.travels(idse,ire)+et0);
        ttmax=ttmin+seconds(mcm.tswind);
        vmin=min(exwfm(:,ire));
        vmax=max(exwfm(:,ire));
        fill([ttmin ttmax ttmax ttmin],[vmin vmin vmax vmax],'r','linestyle','none','FaceAlpha',0.35); hold on;
    end
    xlabel('Time');ylabel('Amplitude');title(sprintf('Station: %s (MCM)',trace.name{ire}));axis tight;
end
seisrsdisp(exwfm,trace.dt,trace.name,exwfm_t0); % display the waveforms of different stations all together

% record section for MCM, without origin time calibration
soup_mcm=search.soup(idse,:)/1000;
if mcm.phasetp==0
    % only P-phase
    dispwfscn(trace.data',trace.recp/1000,soup_mcm,trace.dt,et0,trace.travelp(idse,:),[],[],trace.t0); % note unit transfer
elseif mcm.phasetp==1
    % only S-phase
    dispwfscn(trace.data',trace.recp/1000,soup_mcm,trace.dt,et0,[],trace.travels(idse,:),[],trace.t0); % note unit transfer
elseif mcm.phasetp==2
    % both P- and S-phase
    dispwfscn(trace.data',trace.recp/1000,soup_mcm,trace.dt,et0,trace.travelp(idse,:),trace.travels(idse,:),[],trace.t0); % note unit transfer
end
title('Record section (MCM)');



% Display results, for MCM with origin time calibration
if isfield(mcm,'pperiod')
    tphase=mcm.pperiod;
else
    tphase=0.45; % default period of the seismic phase, in second
end

et0=mcm.st0(tn)+mcm.tpwind-tphase; % the calibrated origin time

% make sure the set time range is not out of boundary
if lwin>et0
    lwinc=et0;
else
    lwinc=lwin;
end
if rwin>mcm.tdatal-et0
    rwinc=mcm.tdatal-et0;
else
    rwinc=rwin;
end

net0r=round((et0-lwinc)/trace.dt+1):round((et0+rwinc)/trace.dt+1); % sample index for extracting waveforms from continuous data
exwfm_t0 = trace.t0 + seconds(et0-lwinc); % t0 time of the extracted data, in datetime format (absolute times)
exwfm=transpose(trace.data(:,net0r)); % extracted waveforms
extmx = trace.t0 + seconds((net0r-1)*trace.dt);  % time axis for waveforms (absolute times)
for ire=1:nre
    figure; plot(extmx,exwfm(:,ire),'k'); hold on;
    if mcm.phasetp==0 || mcm.phasetp==2
        % show P-phase time window
        ttmin=trace.t0+seconds(trace.travelp(idse,ire)+et0);
        ttmax=ttmin+seconds(mcm.tpwind);
        vmin=min(exwfm(:,ire));
        vmax=max(exwfm(:,ire));
        fill([ttmin ttmax ttmax ttmin],[vmin vmin vmax vmax],'b','linestyle','none','FaceAlpha',0.35); hold on;
    end
    if mcm.phasetp==1 || mcm.phasetp==2
        % show S-phase time window
        ttmin=trace.t0+seconds(trace.travels(idse,ire)+et0);
        ttmax=ttmin+seconds(mcm.tswind);
        vmin=min(exwfm(:,ire));
        vmax=max(exwfm(:,ire));
        fill([ttmin ttmax ttmax ttmin],[vmin vmin vmax vmax],'r','linestyle','none','FaceAlpha',0.35); hold on;
    end
    xlabel('Time');ylabel('Amplitude');title(sprintf('Station: %s (MCM t0 calibrated)',trace.name{ire}));axis tight;
end
seisrsdisp(exwfm,trace.dt,trace.name,exwfm_t0); % display the waveforms of different stations all together

% record section for MCM, with origin time calibration
soup_mcm=search.soup(idse,:)/1000;
if mcm.phasetp==0
    % only P-phase
    dispwfscn(trace.data',trace.recp/1000,soup_mcm,trace.dt,et0,trace.travelp(idse,:),[],[],trace.t0); % note unit transfer, m->km
elseif mcm.phasetp==1
    % only S-phase
    dispwfscn(trace.data',trace.recp/1000,soup_mcm,trace.dt,et0,[],trace.travels(idse,:),[],trace.t0); % note unit transfer, m->km
elseif mcm.phasetp==2
    % both P- and S-phase
    dispwfscn(trace.data',trace.recp/1000,soup_mcm,trace.dt,et0,trace.travelp(idse,:),trace.travels(idse,:),[],trace.t0); % note unit transfer, m->km
end
title('Record section (MCM with t0 calibrated)');




if ~isempty(earthquake)
    % Display results, for the catalog
    % determine catalog event index in the soup, by using minimal distance--approximate
    [~,idseca]=min(sum((search.soup-soup_cata).^2,2));
    
    % display the arrival times of seismic event on the recorded seismic data
    et0ca=earthquake.t0;
    net0r=round((et0ca-lwin)/trace.dt+1):round((et0ca+rwin)/trace.dt+1); % origin time for the catalogue
    exwfm_t0 = trace.t0 + seconds(et0-lwinc); % t0 time of the extracted data, in datetime format (absolute times)
    exwfmca=transpose(trace.data(:,net0r)); % extracted waveforms
    for ire=1:nre
        figure; plot((net0r-1)*trace.dt,exwfmca(:,ire),'k'); hold on;
        plot([trace.travelp(idseca,ire)+et0ca   trace.travelp(idseca,ire)+et0ca],[min(exwfmca(:,ire)) max(exwfmca(:,ire))],'b','linewidth',1.2); hold on;
        plot([trace.travels(idseca,ire)+et0ca   trace.travels(idseca,ire)+et0ca],[min(exwfmca(:,ire)) max(exwfmca(:,ire))],'r','linewidth',1.2); hold on;
        xlabel('Time');ylabel('Amplitude');title('Catalogue');axis tight;
    end
    seisrsdisp(exwfmca,trace.dt,trace.name,exwfm_t0); % display the waveforms of different stations all together
    
    % record section of catalogue
    dispwfscn(trace.data',trace.recp/1000,soup_cata/1000,trace.dt,et0ca,trace.travelp(idseca,:),trace.travels(idseca,:),[],trace.t0); % note unit transfer, m->km
    title('Record section (Catalog)');
end


end