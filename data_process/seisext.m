function [data, tt] = seisext(data,dt,t0,timerg)
% This function is used to extract date in the input time range;
% Input:----------------------------------------------------
% data: the recorded seismic data, 2D matrix, dimension: nt*nrec;
% dt: time sampling intveral of the recorded data, unit: s;
% t0: datetime, origin time of the input seismic data;
% timerg: 2*1, time range for extracting date, in datetime format;
%
% OUTPUT-------------------------------------------------------------------
% data: the seismic data in the required time range;


% select date in the defined time range
tid1 = round(seconds(timerg(1) - t0) / dt) + 1;
tid2 = round(seconds(timerg(2) - t0) / dt) + 1;
data = data(tid1:tid2,:);
tt = t0 + seconds(((tid1-1):(tid2-1))*dt);  % time axis for the extracted data

end