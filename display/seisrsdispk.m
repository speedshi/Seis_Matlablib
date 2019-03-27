function seisrsdispk(data,dt)
% This function is used to display the seismic sections and keep the relative amplitude of different traces;
% Display as the record section format according different traces;
% Maximum absolute value of the whole data is linearly normalized to 1,
% so the relative amplitude of each trace is kept.
% Input:----------------------------------------------------
% data: the recorded seismic data, 2D matrix, dimension: nt*nrec
% dt: time intveral for the recorded data, unit: s.

[nt,nrec]=size(data);
drm=bsxfun(@minus,data,mean(data)); % remove the mean values for each trace, i.e. remove the trend of each trace
vmax=max(max(abs(drm)));
if vmax~=0
    dd=drm/vmax;
end
figure;
for ii=1:nrec
    if nargin>1
        plot((0:nt-1)*dt,dd(:,ii)+ii,'k','linewidth',1.1); hold on;
    else
        plot(dd(:,ii)+ii,'k','linewidth',1.1); hold on;
    end
end
if nargin>1
    xlabel('Time (s)');
else
    xlabel('Time points');
end
ylabel('Station number'); axis tight;

end