function migv=wave_migration_kernel_r(trace,mcm,search)
% This function is the MCM kernal.
% The program only saves the maximum migration value at a particular search
% origin time.
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
% mcm.migtp: specify the migration method: 0 for MCM; 1 for DSI;
% search: matlab structure, describe the imaging area;
% search.soup: source imaging positions, 2D array, nsr*3, in meter;
%
% OUTPUT-------------------------------------------------------------------
% migv: migration volume, 2D array, shape: nst0*5; for each row:
% Origin_time-North-East-Depth-Migration_value;


migtp=mcm.migtp; % the migration method: 0 for MCM; 1 for DSI
if migtp==0
    fprintf('Use MCM to locate the earthquake.\n');
elseif migtp==1
    fprintf('Use conventional DSI to locate the earthquake.\n');
else
    error('Incorrect input for mcm.migtp.\n');
end

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

migv=zeros(nst0,5); % initial the final output migration volume

migv_s=zeros(nsr,1); % initial migration volume for one origin time

rsfold='./results'; % name of the output folder
mkdir(rsfold); cd(rsfold); % make and enter the output folder
fname='event_location.dat'; % file name of the ouput file
fid=fopen(fname,'wt'); % open the file for output
fprintf(fid,'%d\n',nst0); % write the total number of origin times into the file


for it=1:nst0
    parfor id=1:nsr
        
        tvpn=round((travelp(id,:)+st0(it))/dt)+1; % time point of direct P-wave for this source position
        tvsn=round((travels(id,:)+st0(it))/dt)+1; % time point of direct S-wave for this source position
        cova_p=zeros(npwd,nre); % initialise the extracted waveform matrix of P-wave, for Z component
        cova_s=zeros(nswd,nre); % initialise the extracted waveform matrix of S-wave, for Z component
        
        for ir=1:nre
            cova_p(:,ir)=data(tvpn(ir):(tvpn(ir)+npwd-1),ir);
            cova_s(:,ir)=data(tvsn(ir):(tvsn(ir)+nswd-1),ir);
        end
        
        if migtp==0
            % use MCM
            pcc=stkcorrcoef(cova_p); % stacked correlation coefficient of P phase
            scc=stkcorrcoef(cova_s); % stacked correlation coefficient of S phase
        else
            % use DSI
            pcc=stkcharfunc(cova_p); % stacked characteristic function of P phase
            scc=stkcharfunc(cova_s); % stacked characteristic function of S phase
        end
        
        migv_s(id)=0.5*(pcc+scc); % the migration value
        
    end
    
    [vmax,indx]=max(migv_s); % obtain the value and index of the maximum migration value
    migv(it,1)=st0(it); % origin time, in second relative to data t0
    migv(it,2:4)=search.soup(indx,:); % potential source locations
    migv(it,5)=vmax; % migration value at this origin time
    
    % write the results to the file
    fprintf(fid,'  %18.6f  %18.6f  %18.6f  %18.6f  %16.8e\n',migv(it,:));
end


fclose(fid);
cd('..');

end