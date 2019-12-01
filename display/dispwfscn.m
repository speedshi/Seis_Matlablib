function dispwfscn(data,recp,soup,dt,stos,travelp,travels,pap,data_t0)
% This function is used to plot the waveform section according to the
% offsets between the event and stations.
% Travel times are marked on the baseline of each traces.
%
% If the input traveltime table is empty, then do not plot it.
%
% Input:----------------------------------------------------
% data: the waveform data recorded at different satations, first dimension:
% time; second dimension: satation.
% recp: the position of different stations, first dimension: station;
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
% data_t0: starting time of seismic data, in matlab datetime format.


% set the default parameters
if ~exist('data_t0','var')
    data_t0=[];
end

[~,nre]=size(data); % obtain the number of time points and stations

% calculate the distances between the event and stations
offset=sqrt(sum(bsxfun(@minus,recp(:,1:3),soup(1,1:3)).^2,2));
[ofsd,idf]=sort(offset); % sort the distance in ascending order, idf is the index for the order
ofdf=diff(ofsd); % offset differences
if min(ofdf)>0.03
    ofsda=ofsd/min(ofdf); % amplify the distances to make the minimal distance interval to be 1
else
    ofsda=ofsd/min(ofdf(ofdf>0.03)); % amplify the distances to make the minimal distance interval to be 1
end

if ~exist('pap','var') || isempty(pap)
    pap=(max(ofsda)-min(ofsda))/nre; % used to amplify the waveform
end

if isempty(travelp) && isempty(travels)
    % both P- and S-traveltimes are empty
    mint=0;
    maxt=0;
elseif ~isempty(travelp) && isempty(travels)
    % have P-traveltimes, no S-traveltimes
    mint=min(travelp); % find the minimal traveltime
    maxt=max(travelp); % find the maximal traveltime
elseif isempty(travelp) && ~isempty(travels)
    % no P-traveltimes, have S-traveltimes
    mint=min(travels); % find the minimal traveltime
    maxt=max(travels); % find the maximal traveltime
elseif ~isempty(travelp) && ~isempty(travels)
    % both P- and S-traveltimes
    mint=min(travelp); % find the minimal traveltime
    maxt=max(travels); % find the maximal traveltime
end

nt1=round((stos+mint-10)/dt)+1; nt2=round((stos+maxt+30)/dt)+1; % obtain the suitable time range for plot the waveform data
% check and make sure the range is not out of boundary
if nt1<1
    nt1=1;
end
if nt2>size(data,1)
    nt2=size(data,1);
end

% plot the waveform section
figure;
for ire=1:nre
    dext=data(nt1:nt2,idf(ire))/max(abs(data(nt1:nt2,idf(ire))))*pap+ofsda(ire); % normalize the waveform data and move the data verticaly according to the distances
    ttx=((nt1-1):(nt2-1))*dt; % time axis
    if isempty(data_t0)
        % no starting time input
        plot(ttx,dext,'k'); hold on; % plot the waveform
    else
        % have starting time input
        plot(data_t0+seconds(ttx),dext,'k'); hold on; % plot the waveform
    end
    if ~isempty(travelp)
        if isempty(data_t0)
            % no starting time input
            plot((travelp(idf(ire))+stos),ofsda(ire),'bx','linewidth',1.2); hold on; % mark the calculated P-wave arrival time
        else
            % have starting time input
            ctime=data_t0+seconds(travelp(idf(ire))+stos);
            plot(ctime,ofsda(ire),'bx','linewidth',1.2); hold on; % mark the calculated P-wave arrival time
        end
    end
    if ~isempty(travels)
        if isempty(data_t0)
            % no starting time input
            plot((travels(idf(ire))+stos),ofsda(ire),'rx','linewidth',1.2); hold on; % mark the calculated S-wave arrival time
        else
            % have starting time input
            ctime=data_t0+seconds(travels(idf(ire))+stos);
            plot(ctime,ofsda(ire),'rx','linewidth',1.2); hold on; % mark the calculated S-wave arrival time
        end
    end
end
axis tight; set(gca,'ytick',ofsda,'yticklabel',num2str(ofsd,'%.2f')); % display the correct distance label
if isempty(data_t0)
    % no starting time input
    xlabel('Time (s)'); ylabel('Distance (km)');
else
    % have starting time input
    xlabel('Time'); ylabel('Distance (km)');
end