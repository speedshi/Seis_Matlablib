function migv=wave_migration_kernel(trace,mcm,search)
% This function is the MCM kernal.
%
% INPUT--------------------------------------------------------------------
% trace: matlab structure, contains seismic data information;
% trace.data: seismic data, 2D array, nre*nt;
% trace.dt: time sampling interval, in second;
% trace.travelp: P-wave traveltime table, 2D array, nsr*nre;
% trace.travels: S-wave traveltime table, 2D array, nsr*nre;
% mcm: matlab structure, contains MCM parameters;
% mcm.tpwind: P-phase time window length in second, scalar;
% mcm.tswind: S-phase time window length in second, scalar;
% mcm.st0: searched origin times of MCM, in second (relative to start time
% of input seismic data), vector, nst0*1;
% search: matlab structure, describe the imaging area,
% search.nsnr: number of imaging points in the north direction, scalar;
% search.nser: number of imaging points in the east direction, scalar;
% search.nsdr: number of imaging points in the depth direction, scalar;
%
% OUTPUT-------------------------------------------------------------------
% migv: migration volume, 4D volume, nsnr*nser*nsdr*nst0.


dt=trace.dt; % time sampling interval of the recorded data (second)

travelp=trace.travelp; % P-wave traveltime table in second
travels=trace.travels; % S-wave traveltime table in second

% format the data, note the data matrix format, it should be: nt*nre
data=trace.data';

st0=mcm.st0; % searched origin times (second)

% migration time window length
npwd=round(mcm.tpwind/dt)+1; % P-wave time window in points
nswd=round(mcm.tswind/dt)+1; % S-wave time window in points


% calculate and set some parameters
[nsr,nre]=size(trace.travelp); % obtain number of source imaging points and stations

nst0=max(size(st0)); % number of searched origin time points

migv=zeros(nsr,nst0); % initial migration volume

parfor it=1:nst0
    for id=1:nsr
        
        tvpn=round((travelp(id,:)+st0(it))/dt)+1; % time point of direct P-wave for this source position
        tvsn=round((travels(id,:)+st0(it))/dt)+1; % time point of direct S-wave for this source position
        cova_p=zeros(npwd,nre); % initialise the extracted waveform matrix of P-wave, for Z component
        cova_s=zeros(nswd,nre); % initialise the extracted waveform matrix of S-wave, for Z component
        
        for ir=1:nre
            cova_p(:,ir)=data(tvpn(ir):(tvpn(ir)+npwd-1),ir);
            cova_s(:,ir)=data(tvsn(ir):(tvsn(ir)+nswd-1),ir);
        end
        
        pcc=stkcorrcoef(cova_p); % stacked correlation coefficient of P phase
        scc=stkcorrcoef(cova_s); % stacked correlation coefficient of S phase
        migv(id,it)=0.5*(pcc+scc); % the final migration value
        
    end
end

% reformat the migration volume, obtain the 4D matrix
migv=reshape(migv,[search.nsnr search.nser search.nsdr nst0]); % reshape to 4D data volume

end