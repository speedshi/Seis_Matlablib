function [exdata,cxdata]=locreson(locres,dt0,ctime,rspac,thrsd)
% This function is used to find events that last a certain period within
% the pre-set distance. Note the unit of the input parameters should keep
% consistent. Time in second (s), distance in meter (m).
% The input events need to be sorted according to their origin time (in
% ascending order).
% INPUT--------------------------------------------------
% locres: event location results, nevt*5, Origin_time-X-Y-Z-Coherency.
% dt0: time interval of the origin times of the located events, in second.
% ctime: pre-set lasting time of an event, only events which last longer
% than this value can be viewed as a true event.
% rspac: pre-set distribution distance of the clustering events, only
% events distributed within this distance can be viewed as lasting in time.
% thrsd: threshold value for the coherency value, only events above the
% threshold value count.
% OUTPUT-----------------------------------------------
% exdata: identified events which fulfill the requirements, ne*5,
% Origin_time-X-Y-Z-Coherency. The position of the events are chosen as the
% same as the events within that time range which have the maximum coherency value.
% oxdata: identified events which fulfill the requirements, ne*5,
% Origin_time-X-Y-Z-Coherency. The position of the events are chosen as the
% central position of the clustering events within that time range.

if nargin<5
    thrsd=0;
end

locres=locres(locres(:,5)>=thrsd,:); % select the events which have coherency value above the threshold

nevt=size(locres,1); % number of seismic events of the input data
nct=round(ctime/dt0); % the limit of lasting time in points

ne=0; % identified events, initialize

% initialize
temp=locres(1,2:4); % used to represent the central position of the clustering events
ncount=1; % the number of lasting time points of events within the distange limit

for ii=2:nevt
    dist=sqrt(sum((locres(ii,2:4)-temp).^2)); % distance between this event and the central position of clustering events
    if dist<rspac && abs(locres(ii,1)-locres(ii-1,1))<=dt0+1e-5
        % within the distance limit, belong to the clustering and can be
        % viewed as the same event
        ncount=ncount+1;
        temp=mean(locres(ii-ncount+1:ii,2:4));
    else
        if ncount>=nct
            % fulfill the lasting time requirement, save this useful event
            ne=ne+1;
            ctemp=locres(ii-ncount:ii-1,:); % pick out the events within that time range
            [~,idx]=max(ctemp(:,5)); % find out the event with maximum coherency value
            exdata(ne,1:5)=ctemp(idx,:);
            exdata(ne,6)=ncount;
            cxdata(ne,:)=exdata(ne,:);
            cxdata(ne,2:4)=temp; % set central position of the clustering events as the output event position
        end
        % move to the next point
        temp=locres(ii,2:4);
        ncount=1;
    end
    
end

end