function [fig, data] = seisrsdisp(data,dt,name,t0,tp,ts,timerg,pp,lw)
% This function is used to display the seismic sections;
% Display as the record section format according different traces;
% Maximum absolute value of each trace is linearly normalized to 1.
% Input:----------------------------------------------------
% data: the recorded seismic data, 2D matrix, dimension: nt*nrec;
% dt: time intveral for the recorded data, unit: s;
% name: station code (name), cell array, shape: nre*1;
% t0: datetime, origin time of the input seismic data;
% tp: nrec*1, P-wave arrivaltime, in datatime format, t0 must be given;
% ts: nrec*1, S-wave arrivaltime, in datetime format, t0 muct be given;
% timerg: 2*1, time range for displaying date, in datetime format, t0 muct be given;
% pp: scalar, for scaling the waveform plots;
% lw: line width for plotting the waveforms;
%
% OUTPUT-------------------------------------------------------------------
% fig: the figure handle;
% data: the seismic data in the required time range;


nrec=size(data,2);

% set default station code and data origin time if no input
if nargin == 2
    name=[];
    t0=[];
    tp=[];
    ts=[];
    timerg=[];
    pp=1;
    lw=1.0;
elseif nargin == 3
    t0=[];
    tp=[];
    ts=[];
    timerg=[];
    pp=1;
    lw=1.0;
elseif nargin == 4
    tp=[];
    ts=[];
    timerg=[];
    pp=1;
    lw=1.0;
elseif nargin == 5
    ts=[];
    timerg=[];
    pp=1;
    lw=1.0;
elseif nargin == 6
    timerg=[];
    pp=1;
    lw=1.0;
elseif nargin == 7
    pp=1;
    lw=1.0;
elseif nargin == 8
    lw=1.0;
end

% if isempty(name)
%     name=1:nrec;
% end

if ~isempty(timerg)
    % select date in the defined time range
    tid1 = round(seconds(timerg(1) - t0) / dt) + 1;
    tid2 = round(seconds(timerg(2) - t0) / dt) + 1;
    data = data(tid1:tid2,:);
    t0 = t0 + seconds((tid1-1)*dt);
end

nt = size(data,1);
fig=figure;
for ii=1:nrec
    
    drm=data(:,ii)-mean(data(:,ii)); % remove the mean values for each trace, i.e. remove the trend of each trace
    
    % shift the data
    if max(abs(drm))>eps
        dd=drm/max(abs(drm))*pp+ii;
    else
        dd=drm+ii;
    end
    
    if nargin==1
        plot(dd,'k','linewidth',lw); hold on;
    elseif nargin<4
        plot((0:nt-1)*dt,dd,'k','linewidth',lw); hold on;
    else
        xtime=t0+seconds((0:nt-1)*dt);
        plot(xtime,dd,'k','linewidth',lw); hold on;
    end
    
    % plot P-wave arrivaltime
    if ~isempty(tp)
        plot(tp(ii),ii,'bx','linewidth',1.1);
    end
    
    % plot S-wave arrivaltime
    if ~isempty(ts)
        plot(ts(ii),ii,'rx','linewidth',1.1);
    end
    
end

% set xlabel
if nargin==1
    xlabel('Time points');
elseif nargin<4
    xlabel('Time (s)');
else
    xlabel('Time');
end

% set ylable
if ~isempty(name)
    set(gca,'ytick',1:nrec,'yticklabel',name);
end
ylabel('Station code');axis tight;

end
