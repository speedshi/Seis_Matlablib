function fdata=my_kurtosis(data,twn,dtp)
% This function is used to calculate the Kurtosis of the seismogram within
% a specified time window.
% Input:----------------------------------------------------
% 1 data: input seismogram, could be a trace (1 dimension) or a recorded
% profile (2 dimension). The seismic time sampling data must be stored in
% the 1st dimension;
% 2 twn: length of time window in samples;
% 3 dtp: set the type of input data to calculate Kurtosis (0 for using the
% square value to calculate Kurtosis; 1 for using absolute value to
% calculate Kurtosis; 2 for using envolope to calculate Kurtosis; 3 for using
% original data to calculate Kurtosis-default).
% Output:--------------------------------------------------
% fdata: the calculated kurtosis of the seismogram within the selected time
% window; has the same dimension with input 'data'.

if nargin<3
    dtp=3;
end

NT=size(data,1); % the number of time samples

% transfer the input data to particular type
if dtp==0
    % square value
    tdata=data.^2;
elseif dtp==1
    % absolute value
    tdata=abs(data);
elseif dtp==2
    % envolope
    tdata=abs(hilbert(data));
elseif dtp==3
    % original data
    tdata=data;
else
    error('Wrong input for dtp!');
end

fdata=zeros(size(data)); % initialize the output

for it=twn:NT
    fdata(it,:)=kurtosis(tdata(it-twn+1:it,:));
end

end