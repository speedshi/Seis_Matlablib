function sdata=wave_extract(data,dt,atimes,twinl,twinr)
% This function is used to extract waveforms from the input data according
% to the input arrival times. This is used to allign the traces for different
% stations.
%
% The arrival times 'atimes' must be expressed in second relative to the
% starting time of the traces.
%
% INPUT--------------------------------------------------------------------
% data: input seismic data, 2D array, shape: nt*nre;
% dt: time sample interval, scalar, in second;
% atimes: arrival times for differen stations in second, 1D array, shape: nre*1;
% twinl: left time window for extract the waveforms, in second, scalar;
% twinr: right time window for extract the waveforms, in second, scalar;
%
% OUTPUT-------------------------------------------------------------------
% sdata: extracted seismic data, 2D array, shape: nt_e*nre;

[~,nre]=size(data); % obtain the number of time points and stations

tpl=round(twinl/dt); % time points on the left
tpr=round(twinr/dt); % time points on the right

sdata=zeros(tpl+tpr+1,nre); % initial the output array

for ii=1:nre
    tindx=round(atimes(ii)/dt)+1; % the sample index of the arrival time
    sdata(:,ii)=data(tindx-tpl:tindx+tpr,ii); % extract the data
end



end