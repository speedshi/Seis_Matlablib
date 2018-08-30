function fdata=intder(data,dt,mtp)
% This function is used to calculate the integration or derivative of the
% input singal.
% Input:----------------------------------------------------
% 1 data: input seismogram (signal), could be a trace (1 dimension) or a recorded
% profile (2 dimension). The seismic time sampling data must be stored in
% the 1st dimension;
% 2 dt: time interval of the data (s);
% 3 mtp: 0 for calculating the derivative of the input data; 1 for
% calculating the integration of the input data.
% Output:--------------------------------------------------
% fdata: the calculated integration or derivative of the input data; has
% the same size with the input 'data'.
% Note we use central difference to calculate the derivatives and integrations.
% For integration, integration starts from the first time point.


NT=size(data,1); % number of time points

fdata=zeros(size(data)); % initialize the output data
if mtp==0
    % calculate the derivatives using central difference
    for it=2:NT-1        
        fdata(it,:)=0.5*(data(it+1,:)-data(it-1,:))/dt;
    end
elseif mtp==1
    % calculate the integrations using central difference
    fdata(1,:)=data(1,:)*dt;
    for it=2:NT
        fdata(it,:)=fdata(it-1,:)+(data(it,:)+data(it-1,:))*dt/2;
    end   
end

end