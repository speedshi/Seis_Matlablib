function fstack=dispwflstk(data,dt,stos,travelt,tw1,tw2)
% This function is used to plot the linear waveform stack results according
% to the arrival times. The figure is used to check whether there is any
% signals at the input arrival times.
% INPUT-------------------------------------------------------------------------
% data: the waveform data recorded at different satations, nt*nre, first
% dimension: time; second dimension: satation.
% dt: time sampling interval of the recorded data, in second.
% travelt: the travel times of a particular phase from this event to
% different stations, 1*nre, different column corresponding to different
% stations, in second.
% tw1: the left time window length used for stacking, in second.
% tw2: the right time window length used for stacking, in second.

if nargin<=4
    tw1=1; % the left time window used for stacking
    tw2=1; % the right time window used for stacking
end

nre=length(travelt); % number of stations

ntl=(tw2+tw1)/dt+1;
fstack=zeros(ntl,1);
stn=tw1/dt+1;

for ii=1:nre
    % obtain the suitable time range for stacking waveform data
    nt1=round((stos+travelt(ii)-tw1)/dt)+1;
    nt2=round((stos+travelt(ii)+tw2)/dt)+1;
    snn=round((stos+travelt(ii))/dt)+1;
    if sum(data(snn:snn+4,ii))>0
       scoe=-1;
    else
        scoe=1; % station weighting coefficient, used to account for source radiation pattern
    end
    fstack=fstack+scoe*data(nt1:nt2,ii)/max(abs(data(nt1:nt2,ii)));
end

fstack=fstack/nre;

figure; plot(-tw1:dt:tw2,fstack,'k'); hold on; % plot the waveform
plot(0,fstack(stn),'rx');
axis tight; xlabel('Time (second)'); ylabel('Amplitude');

end