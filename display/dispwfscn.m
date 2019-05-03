function dispwfscn(data,recp,soup,dt,stos,travelp,travels,pap)
% This function is used to plot the waveform section according to the
% offsets between the event and stations.
% Travel times are marked on the baseline of each traces.
% Input:----------------------------------------------------
% data: the waveform data recorded at different satations, first dimension:
% time; second dimension: satation.
% recp: the position of different stations, first dimension: station;
% second dimension: position (column 1-3: X-Y-Z), in km.
% soup: the position of a particular event, 1*3, column 1-3: X-Y-Z, in km.
% dt: time sampling interval of the recorded data, in second.
% stos: the original time of this event relative to the initial time of the
% recorded data, in second.
% travelp: the travel time of P-wave from this event to different stations,
% 1*nre, different column corresponding to different stations, in second.
% travels: the travel time of S-wave from this event to different stations,
% 1*nre, different column corresponding to different stations, in second.
% pap: used to amplify the waveform.

[~,nre]=size(data); % obtain the number of time points and stations

% calculate the distances between the event and stations
offset=sqrt(sum(bsxfun(@minus,recp(:,1:3),soup(1,1:3)).^2,2));
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

mint=min(travelp); % find the minimal travel time
maxt=max(travels); % find the maximal travel time

nt1=round((stos+mint-5)/dt)+1; nt2=round((stos+maxt+8)/dt)+1; % obtain the suitable time range for plot the waveform data

% plot the waveform section
figure;
for ire=1:nre
    dext=data(nt1:nt2,idf(ire))/max(abs(data(nt1:nt2,idf(ire))))*pap+ofsda(ire); % normalize the waveform data and move the data verticaly according to the distances
    plot(((nt1-1):(nt2-1))*dt,dext,'k'); hold on; % plot the waveform
    plot((travelp(idf(ire))+stos),ofsda(ire),'bx','linewidth',1.2); hold on; % mark the calculated P-wave arrival time
    plot((travels(idf(ire))+stos),ofsda(ire),'rx','linewidth',1.2); hold on; % mark the calculated S-wave arrival time
end
axis tight; set(gca,'ytick',ofsda,'yticklabel',num2str(ofsd)); % display the correct distance label
xlabel('Time (s)'); ylabel('Distance (km)');