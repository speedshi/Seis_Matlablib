function migv=wave_migration_kernel_x(trace,mcm,search)
% This function is the MCM kernal.
% Note this function only use single phase.
%
% INPUT--------------------------------------------------------------------
% trace: matlab structure, contains seismic data information;
% trace.data: seismic data, 2D array, nre*nt;
% trace.dt: time sampling interval, in second;
% trace.travelx: traveltime table of the seismic phase, 2D array, nsr*nre;
% mcm: matlab structure, contains MCM parameters;
% mcm.txwind: time window length in second for the seismic phase, scalar;
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

travelx=trace.travelx; % traveltime table in second of the seismic phase

% format the data matrix, note the data matrix format, it should be: nt*nre
% in this function
data=trace.data';

st0=mcm.st0; % searched origin times (second)

% migration time window length
nxwd=round(mcm.txwind/dt)+1; % time window in points for the seismic phase


% calculate and set some parameters
[nsr,nre]=size(trace.travelx); % obtain number of source imaging points and stations

ncoe=nre*(nre-1)/2; % number of unique station pairs

nst0=max(size(st0)); % number of searched origin time points

migv=zeros(nsr,nst0); % initial migration volume

parfor it=1:nst0
    for id=1:nsr
        
        tvxn=round((travelx(id,:)+st0(it))/dt)+1; % time point of the seismic phase for this source position
        cova_x=zeros(nxwd,nre); % initialise the extracted waveform matrix of the seismic phase
        
        for ir=1:nre
            cova_x(:,ir)=data(tvxn(ir):(tvxn(ir)+nxwd-1),ir);
        end
        
        xcc=triu(mycorrcoef(cova_x),1); % correlation coefficient matrix of the seismic phase
        
        migv(id,it)=sum(abs(xcc(:)))/(ncoe);
        
    end
end

% reformat the migration volume, obtain the 4D matrix
migv=reshape(migv,[search.nsnr search.nser search.nsdr nst0]); % reshape to 4D data volume

end