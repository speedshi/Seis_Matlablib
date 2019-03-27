function dispwfscn_2se(data,recp,soup,dt,stos,travelp,travels,pap)
% This function is used to plot the waveform section according to the
% offsets between the event and stations. Here we input 2 source events
% (the same event located by different methods, this record section is used
% to compare the differnces of these two located events) and display the
% predicted arrivals of the two events at the record section. The offsets
% are the average offsets for event1 offsets and event2 offsets.
% Input:----------------------------------------------------
% data: the waveform data recorded at different satations, first dimension:
% time; second dimension: satation.
% recp: the position of different stations, first dimension: station;
% second dimension: position (column 1-3: X-Y-Z), in km.
% soup: the position of a particular event, 2*3, column 1-3: X-Y-Z, in km;
% row 1: event 1, row 2: event 2.
% dt: time sampling interval of the recorded data, in second.
% stos: the original time of the 2 event relative to the initial time of the
% recorded data, in second, 2*1 vector.
% travelp: the travel time of P-wave from the 2 events to different stations,
% 2*nre, different column corresponding to different stations, in second.
% travels: the travel time of S-wave from the 2 events to different stations,
% 2*nre, different column corresponding to different stations, in second.
% pap: used to amplify the waveform.

[~,nre]=size(data); % obtain the number of time points and stations

% calculate the distances between the event and stations
offset1=sqrt(sum(bsxfun(@minus,recp(:,1:3),soup(1,1:3)).^2,2)); % offsets of event 1
offset2=sqrt(sum(bsxfun(@minus,recp(:,1:3),soup(2,1:3)).^2,2)); % offsets of event 2
offset=0.5*(offset1+offset2); % average offsets
[ofsd,idf]=sort(offset); % sort the distance in ascending order, idf is the index for the order
ofdf=diff(ofsd); % offset differences
if min(ofdf)>0.03
    ofsda=ofsd/min(ofdf); % amplify the distances to make the minimal distance interval to be 1
else
    ofsda=ofsd/min(ofdf(ofdf>0.03)); % amplify the distances to make the minimal distance interval to be 1
end

if nargin<8
    pap=(max(ofsda)-min(ofsda))/nre; % used to amplify the waveform
end

mint1=min(travelp(1,:)); % find the minimal travel time of event 1
maxt1=max(travels(1,:)); % find the maximal travel time of event 1
mint2=min(travelp(2,:)); % find the minimal travel time of event 2
maxt2=max(travels(2,:)); % find the maximal travel time of event 2

nt1=round(min(stos(1)+mint1-1,stos(2)+mint2-1)/dt)+1;
nt2=round(max(stos(1)+maxt1+3,stos(2)+maxt2+3)/dt)+1; % obtain the suitable time range for plot the waveform data

% plot the waveform section
figure;
plot((travelp(1,idf)+stos(1)),ofsda,'-b','linewidth',1.2); hold on; % mark the calculated P-wave arrival time for event 1
plot((travelp(2,idf)+stos(2)),ofsda,'--b','linewidth',1.2); hold on; % mark the calculated P-wave arrival time for event 2
for ire=1:nre
    dext=data(nt1:nt2,idf(ire))/max(abs(data(nt1:nt2,idf(ire))))*pap+ofsda(ire); % normalize the waveform data and move the data verticaly according to the distances
    plot(((nt1-1):(nt2-1))*dt,dext,'k'); hold on; % plot the waveform
    plot((travelp(1,idf(ire))+stos(1)),ofsda(ire),'bx','linewidth',1.2); hold on; % mark the calculated P-wave arrival time for event 1
    plot((travels(1,idf(ire))+stos(1)),ofsda(ire),'rx','linewidth',1.2); hold on; % mark the calculated S-wave arrival time for event 1
    plot((travelp(2,idf(ire))+stos(2)),ofsda(ire),'bx','linewidth',1.2); hold on; % mark the calculated P-wave arrival time for event 2
    plot((travels(2,idf(ire))+stos(2)),ofsda(ire),'rx','linewidth',1.2); hold on; % mark the calculated S-wave arrival time for event 2
end
plot((travels(1,idf)+stos(1)),ofsda,'-r','linewidth',1.2); hold on; % mark the calculated S-wave arrival time for event 1
plot((travels(2,idf)+stos(2)),ofsda,'--r','linewidth',1.2); hold on; % mark the calculated S-wave arrival time for event 2
axis tight; set(gca,'ytick',ofsda,'yticklabel',num2str(ofsd)); % display the correct distance label
xlabel('Time (s)'); ylabel('Distance (km)');