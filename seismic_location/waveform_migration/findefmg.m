function events=findefmg(data,thrsd,eots)
% This function is used to find the real events from the preliminary migration results.
% The input events need to be sorted according to their origin time (in
% ascending order). The unit of origin time should be second.
% INPUT--------------------------------------------------
% data: event location results, nevt*5, Origin_time-X-Y-Z-Coherency.
% thrsd: threshold value, only events above the threshold are viewed as real events.
% eots: event origin time separation limit, in second. Within this time limit,
% they are viewed as the same event; larger than this limit, they are
% viewed as different events.
% OUTPUT-----------------------------------------------
% events: the final found real events, ne*5, Origin_time-X-Y-Z-Coherency.

sdata=data(data(:,5)>=thrsd,:); % find events which have coherency value larger than the threshold
npet=size(sdata,1); % the number of selected potential events

if npet<2
    % no or only one event found
    events=sdata;
    return;
end

ne=1; % number of the found real events, initialize
events(ne,:)=sdata(1,:); % initialize the output results

for ii=2:npet
    if sdata(ii,1)-events(ne,1)<=eots % sdata(ii,1)-sdata(ii-1,1)<=eots
        % within the origin time separation limit, belong to the same event
        if sdata(ii,5)>events(ne,5)
            % replace the current event using the event which has higher coherency
            events(ne,:)=sdata(ii,:);
        end
    else
        % out of the separation limit, move to another event
        ne=ne+1;
        events(ne,:)=sdata(ii,:);
    end
end

end