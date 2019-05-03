function data=extractevt(cata,mloc,t1,t2)
% This function is used to extract seismic events within the time range of
% reference seismic events in the catalogue. Identifying the event with
% maximum coherency value within the reference time range.
% Notice the unit of time is second (s).
% INPUT:----------------------------------------------------
% cata: reference seismic event locations in the catalogue, neca*4, Origin_time-X-Y-Z;
% mloc: our event location results, nevt*5, Origin_time-X-Y-Z-Coherency;
% t1: the left accepted time rang (limit) in second (s);
% t2: the right accepted time rang (limit) in second (s).
% OUTPUT-----------------------------------------------
% data: identified seismic events within the referenced time range, neca*5,
% Origin_time-X-Y-Z-Coherency.

% set default value
if nargin<=2
    t1=1;
    t2=0;
end

neca=size(cata,1); % the number of seismic events in the input catalogue
data=zeros(neca,5);

for ii=1:neca
    % set the time range
    tmin=cata(ii,1)-t1; % low boundary of time range
    tmax=cata(ii,1)+t2; % upper boundary of time range
    
    tempa=mloc(mloc(:,1)>=tmin & mloc(:,1)<=tmax,:); % find located seismic events in this time range
    [~,indx]=max(tempa(:,5)); % find the event with maximum coherency value
    data(ii,:)=tempa(indx,:);
end

end