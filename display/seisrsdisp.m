function seisrsdisp(data,dt)
% This function is used to display the seismic sections;
% Display as the record section format according different traces;
% Maximum absolute value of each trace is linearly normalized to 1.
% Input:----------------------------------------------------
% data: the recorded seismic data, 2D matrix, dimension: nt*nrec
% dt: time intveral for the recorded data, unit: s.

[nt,nrec]=size(data);
figure;
for ii=1:nrec
    drm=data(:,ii)-mean(data(:,ii)); % remove the mean values for each trace, i.e. remove the trend of each trace
    if max(abs(drm))>eps
        dd=drm/max(abs(drm))+ii;
    else
        dd=drm+ii;
    end
    if nargin>1
        plot((0:nt-1)*dt,dd,'k','linewidth',1.1); hold on;
    else
        plot(dd,'k','linewidth',1.1); hold on;
    end
end
if nargin>1
    xlabel('Time (s)');
else
    xlabel('Time points');
end
ylabel('Station number'); axis tight;

end