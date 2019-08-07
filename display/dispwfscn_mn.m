function dispwfscn_mn(data,recp,soup,dt,stos,travelp,travels,name,t0,pap)
% This function is used to plot the waveform section according to the
% offsets between the event and stations. The traces are separated and
% sorted according to source-receiver distances and normalized.
% Arrival-times are marked on the baseline of each traces.
% Station names are shown in the figure accordingly.
% Input:----------------------------------------------------
% data: the waveform data, 2D array, shape: nt*nre;
% recp: the position of stations, shape: nre*3, in km; North-East-Depth;
% soup: the position of a particular event, 1*3, North-East-Depth, in km.
% dt: time sampling interval of the recorded data, in second.
% stos: the original time of this event relative to the initial time of the
% recorded data, in second.
% travelp: the travel time of P-wave from this event to different stations,
% 1*nre, different column corresponding to different stations, in second.
% travels: the travel time of S-wave from this event to different stations,
% 1*nre, different column corresponding to different stations, in second.
% pap: used to amplify the waveform;
% name: names of the stations, cell array, shape: nre*1;
% t0: matlab datetime, the origin time of the seismic data;

[~,nre]=size(data); % obtain the number of time points and stations

% calculate the distances between the event and stations
offset=sqrt(sum(bsxfun(@minus,recp(:,1:3),soup(1,1:3)).^2,2));

[ofsd,idf]=sort(offset); % sort the distance in ascending order, idf is the index for the order


mint=min(travelp); % find the minimal travel time
maxt=max(travels); % find the maximal travel time

nt1=round((stos+mint-10)/dt)+1; nt2=round((stos+maxt+30)/dt)+1; % obtain the suitable time range for plot the waveform data

% plot the waveform section
figure;
for ire=1:nre
    dext=data(nt1:nt2,idf(ire))/max(abs(data(nt1:nt2,idf(ire))))+ire; % normalize the waveform data and move the data verticaly accordingly
    xtimes=t0+seconds(((nt1-1):(nt2-1))*dt); % time axis
    plot(xtimes,dext,'k'); hold on; % plot the waveform
    xtimes=t0+seconds(travelp(idf(ire))+stos); % time axis
    plot(xtimes,ire,'bx','linewidth',1.2); hold on; % mark the calculated P-wave arrival time
    xtimes=t0+seconds(travels(idf(ire))+stos); % time axis
    plot(xtimes,ire,'rx','linewidth',1.2); hold on; % mark the calculated S-wave arrival time
end
axis tight; set(gca,'ytick',1:nre,'yticklabel',name(idf)); % display the correct distance label
xlabel('Time'); ylabel('Station code');
