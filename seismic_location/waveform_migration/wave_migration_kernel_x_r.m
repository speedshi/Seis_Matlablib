function migv=wave_migration_kernel_x_r(trace,mcm,search)
% This function is the MCM kernal.
% Note this function only use single phase.
% The program only saves the maximum migration value at a particular search
% origin time.
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

travelx=trace.travelx; % traveltime table in second of the seismic phase

% format the data matrix, note the data matrix format, it should be: nt*nre
% in this function
data=trace.data';

st0=mcm.st0; % searched origin times (second)

% migration time window length
nxwd=round(mcm.txwind/dt)+1; % time window in points for the seismic phase


% calculate and set some parameters
[nsr,nre]=size(trace.travelx); % obtain number of source imaging points and stations

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
        
        tvxn=round((travelx(id,:)+st0(it))/dt)+1; % time point of the seismic phase for this source position
        cova_x=zeros(nxwd,nre); % initialise the extracted waveform matrix of the seismic phase
        
        for ir=1:nre
            cova_x(:,ir)=data(tvxn(ir):(tvxn(ir)+nxwd-1),ir);
        end
        
        if migtp==0
            % use MCM
            migv_s(id)=stkcorrcoef(cova_x); % stacked correlation coefficient of the seismic phase
        else
            % use DSI
            migv_s(id)=stkcharfunc(cova_x); % stacked characteristic function of the seismic phase
        end
        
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