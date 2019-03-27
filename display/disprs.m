function [ofsda,ofsd]=disprs(data,recp,soup,dt,pap,dpre,stos,travelp,travels)
% This function is used to plot the record section according to the
% horizontal offsets between the event and stations.
% Input:----------------------------------------------------
% data: the waveform data recorded at different satations, first dimension:
% time; second dimension: satation.
% recp: the position of different stations, Nre*3, first dimension: station;
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
% Output:--------------------------------------------------
% ofsda: the amplified horizontal offset (km) used for Y-axis.
% ofsd: the true horizontal offsets (km) of the plotted traces.


if nargin<5
    pap=10; % used to amplify the waveform
    dpre=3; % set the precision of distance in km, expressed as the number after decimal point
    iflagt=0; % used to contral the plot time range
    iflagp=0; % used to contral if plot the P-wave arrivals
    iflags=0; % used to contral if plot the S-wave arrivals
elseif nargin<6
    dpre=3; % set the precision of distance in km, expressed as the number after decimal point
    iflagt=0; % used to contral the plot time range
    iflagp=0; % used to contral if plot the P-wave arrivals
    iflags=0; % used to contral if plot the S-wave arrivals
elseif nargin<7
    iflagt=0; % used to contral the plot time range
    iflagp=0; % used to contral if plot the P-wave arrivals
    iflags=0; % used to contral if plot the S-wave arrivals
elseif nargin<9
    iflagt=1; % used to contral the plot time range
    iflagp=1; % used to contral if plot the P-wave arrivals
    iflags=0; % used to contral if plot the S-wave arrivals
else
    iflagt=1; % used to contral the plot time range
    iflagp=1; % used to contral if plot the P-wave arrivals
    iflags=1; % used to contral if plot the S-wave arrivals
end

[Nt,nre]=size(data); % obtain the number of time points and stations

% calculate the distances between the event and stations
offseto=sqrt(sum(bsxfun(@minus,recp(:,1:2),soup(1,1:2)).^2,2));
offset=round(offseto,dpre); % set the precision of the source-receiver offsets,
[ofsd,idf,~]=unique(offset); % sort the distance in ascending order, idf is the index for the order
if min(diff(ofsd))>0
    ofsda=ofsd/min(diff(ofsd)); % amplify the distances to make the minimal distance interval to be 1
else
    ofsda=ofsd;
end

if iflagt==1
    mint=min(travelp); % find the minimal travel time
    maxt=max(travels); % find the maximal travel time
    nt1=round((stos+mint-1)/dt)+1; nt2=round((stos+maxt+3)/dt)+1; % obtain the suitable time range for plot the waveform data
    if nt1<1
        nt1=1;
    end
    if nt2>Nt
        nt2=Nt;
    end
else
    nt1=1; nt2=Nt;
end

% plot the waveform section
figure;
for ire=1:length(ofsd)
    vnorm=max(abs(data(nt1:nt2,idf(ire))));
    if vnorm==0
        vnorm=1;
    end
    dext=data(nt1:nt2,idf(ire))/vnorm*pap+ofsda(ire); % normalize the waveform data and move the data verticaly according to the distances
    plot(((nt1-1):(nt2-1))*dt,dext,'k'); hold on; % plot the waveform
    if iflagp==1
        plot((travelp(idf(ire))+stos),dext(round((travelp(idf(ire))+stos)/dt)+2-nt1),'b.','linewidth',1.2); hold on; % mark the calculated P-wave arrival time
    end
    if iflags==1
        plot((travels(idf(ire))+stos),dext(round((travels(idf(ire))+stos)/dt)+2-nt1),'r.','linewidth',1.2); hold on; % mark the calculated S-wave arrival time
    end
end
plot((travelp(idf)+stos),ofsda,'-b','linewidth',1.2); hold on; % mark the calculated P-wave arrival time - line
plot((travels(idf)+stos),ofsda,'-r','linewidth',1.2); hold on; % mark the calculated S-wave arrival time - line
axis tight; set(gca,'ytick',ofsda,'yticklabel',num2str(ofsd)); % display the correct distance label
xlabel('Time (s)'); ylabel('Offset (km)');