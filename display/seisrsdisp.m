function seisrsdisp(data,dt,name,t0)
% This function is used to display the seismic sections;
% Display as the record section format according different traces;
% Maximum absolute value of each trace is linearly normalized to 1.
% Input:----------------------------------------------------
% data: the recorded seismic data, 2D matrix, dimension: nt*nrec;
% dt: time intveral for the recorded data, unit: s;
% name: station code (name), cell array, shape: nre*1;
% t0: datetime, origin time of the input seismic data;


[nt,nrec]=size(data);

% set default station code and data origin time if no input
if nargin == 2
    name=1:nrec;
    t0=[];
end

figure;
for ii=1:nrec
    drm=data(:,ii)-mean(data(:,ii)); % remove the mean values for each trace, i.e. remove the trend of each trace
    
    % shift the data
    if max(abs(drm))>eps
        dd=drm/max(abs(drm))+ii;
    else
        dd=drm+ii;
    end
    
    if nargin>1 && nargin<4
        plot((0:nt-1)*dt,dd,'k','linewidth',1.1); hold on;
    elseif nargin==4
        xtime=t0+seconds((0:nt-1)*dt);
        plot(xtime,dd,'k','linewidth',1.1); hold on;
    else
        plot(dd,'k','linewidth',1.1); hold on;
    end
end

if nargin>1 && nargin<4
    xlabel('Time (s)');
elseif nargin==4
    xlabel('Time');
else
    xlabel('Time points');
end

set(gca,'ytick',1:nrec,'yticklabel',name);
ylabel('Station code');axis tight;

end