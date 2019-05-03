function gpltlocrs(data,events,time,cs)
% This function is used to plot the corresponding located events.
% Used for evaluation of the location performance.
% INPUT--------------------------------------------------
% data: whole migration data,  nevt*5, Origin_time-X-Y-Z-Coherency;
% events: the found real events, ne*5, Origin_time-X-Y-Z-Coherency;
% time: the time axis of the data, matlab datatime format;
% cs: the color of the highlighted events.

if nargin==2
    time=[];
    cs='r';
elseif nargin==3
    cs='r';
end

ne=size(events,1); % number of real events
sid=zeros(ne,1); % events index

for ii=1:ne
    indx=find(abs(data(:,1)-events(ii,1))<1e-6); % find the events index according its origin time, origin time is unique
    if isempty(indx)
        sid(ii)=NaN;
    else
        sid(ii)=indx;
    end
end

ssid=sid(~isnan(sid)); % remove the NAN values in the index

if isempty(time)
    plot(data(:,1),data(:,5),'k'); hold on;
    plot(data(ssid,1),data(ssid,5),[cs 'o'],'markersize',3,'markerfacecolor',cs); hold on;
    xlim([data(ssid(1),1) data(ssid(end),1)]); xlabel('Time'); ylabel('Coherency');
    %set(gca,'xtick',[0 3600 7200 10800 14400],'xticklabel',{'06:00:00','07:00:00','08:00:00','09:00:00','10:00:00'});
else
    plot(time,data(:,5),'k'); hold on;
    plot(time(ssid),data(ssid,5),[cs 'o'],'markersize',3,'markerfacecolor',cs); hold on;
    xlim([time(ssid(1),1) time(ssid(end),1)]);xlabel('Time'); ylabel('Coherency');
end


end