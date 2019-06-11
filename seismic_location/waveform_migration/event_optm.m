function cdata=event_optm(migloc,timelim,spaclim,vthrd)
% This function is used to optimise the input event location results.
% Unit should keep consistent for the data in: 'cata' and parameters: 'timelim' and 'spaclim'.
% Unit should use second and meter.
%
% Input:----------------------------------------------------
% migloc: the original location results, shape: nevt*5, in Time-North-East-Depth-Coherency;
%
% timelim: in second, minimal time limit, indicates the close located
% seismic events must have a time interval larger than 'timelim';
%
% spaclim: in meter, minimal space limit, indicates the close located
% seismic events must have a space interval larger than 'spaclim';
%
% vthrd: only events with coherency value >= "vthrd" are choosed as seismic events.
%
% Output:--------------------------------------------------
% cdata: event locations after optimization, neo*5, Time-North-East-Depth-Coherency.

% set default value
if nargin<=3
    vthrd=0;
end

data_temp=migloc(migloc(:,5)>=vthrd,:); % choose events which have a coherency value >= pre-set threshold

% sort the events according their coherency values in descending order
[~,indx]=sort(data_temp(:,5),'descend');
data_temp=data_temp(indx,:);

npe=size(data_temp,1); % the total number of potential events
if npe>0
    neo=1; % number of effective seismic events that has been found
    cdata(1,:)=data_temp(1,:); % the event with the largest coherency value must be an event
else
    cdata=[0 0 0 0 0];
    fprintf('Threshold value is too large! No events exceed the threshold!\n');
    return;
end

for ii=2:npe
    iflag=1; % used to show is this an event or not
    for ie=1:neo
        dist=sqrt(sum((data_temp(ii,2:4)-cdata(ie,2:4)).^2)); % space interval between two events
        tnit=abs(data_temp(ii,1)-cdata(ie,1)); % time interval between two events
        if (dist<=spaclim && tnit<=timelim)
            % within the space and time linit, not an effective event
            iflag=0;
            break;
        end
    end
    
    if iflag==1
        % fulfill the space and time interval requirement, is an effective seismic event
        neo=neo+1;
        cdata(neo,:)=data_temp(ii,:);
    end
end

end