function fdata=my_stalta(data,stn,ltn,dtp,mtp)
% This function is used to calculate the STA/LTA ratio.
% Input:-----------------------------------------------------------
% 1 data: seismic data in which the 1st dimension contains the time
% sampling data;
% 2 stn: Length of short time average window in samples;
% 3 ltn: Length of long time average window in samples;
% 4 dtp: set the type of input data to calculate STA/LTA (0 for using the
% square value to calculate STA/LTA; 1 for using absolute value to
% calculate STA/LTA; 2 for using envolope to calculate STA/LTA; 3 for using
% the output of Z-detector to calculate STA/LTA);
% 5 mtp: set the type of method to calculate STA/LTA (0 for calculating the
% classical STA/LTA; 1 for calculate the recursive STA/LTA; 2 for calculate
% the delayed STA/LTA; 3 for calcualting the classical STA/LTA with
% non-overlapping time window; 4 no STA/LTA calculation, output the
% transferred data directly).
% Output:--------------------------------------------------------
% fdata: the calculated STA/LTA value in every time point.
% An exhibition of time window for classical STA/LTA:
% ttttttttttttttttttttttttttttttttttttttttttttttTtttttttttttttttt
%                               |<--STA-->|
% |<---------------LTA-------------->|
% An exhibition of time window for classical STA/LTA with non-overlapping time window:
% ttttttttttttttttttttttttttttttttttttttttttttttTtttttttttttttttt                              
% |<---------------LTA-------------->|<--STA-->|

% Set default value
if nargin<4
    dtp=0;
    mtp=0;
elseif nargin<5
    mtp=0;
end

[NT,nre]=size(data); % Nt: number of time samples; nre: the number of receivers

tdata=zeros(size(data));
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
    % Z-detector
    zsd=zeros(size(data));
    for it=stn:NT
        zsd(it,:)=sum(data(it-stn+1:it,:));
    end
    z_mean=mean(zsd);
    z_std=std(zsd);
    for ir=1:nre
        tdata(:,ir)=(zsd(:,ir)-z_mean(ir))/z_std(ir);
    end
else
    error('Wrong input for dtp!');
end

fdata=zeros(size(data)); % initialize the output data
% the time samples at the begining and end are set to be 0.


if mtp==0
    % calculate the classical STA/LTA
    for it=ltn:NT
        sta_ave=mean(tdata(it-stn+1:it,:)); % short time average
        lta_ave=mean(tdata(it-ltn+1:it,:)); % long time average
        lta_ave(lta_ave<eps)=eps; % avoid divided by 0
        fdata(it,:)=sta_ave./lta_ave;
    end
elseif mtp==1
    % calculate the recursive STA/LTA
    stas=zeros(1,size(data,2));
    ltas=eps*ones(1,size(data,2)); % avoid divided by 0
    csta=1/stn;
    clta=1/ltn;
    for it=2:NT
        stas=csta*tdata(it,:)+(1-csta)*stas;
        ltas=clta*tdata(it,:)+(1-clta)*ltas;
        fdata(it,:)=stas./ltas;
        if it<=ltn
            fdata(it,:)=0; % set 0 within the long time window
        end
    end
elseif mtp==2
    % calculate the delayed STA/LTA
    sta=zeros(size(data));
    lta=zeros(size(data));
    for it=stn+ltn+2:NT
        sta(it,:)=(tdata(it,:)+tdata(it-stn,:))/stn+sta(it-1,:);
        lta(it,:)=(tdata(it-stn-1,:)+tdata(it-stn-ltn-1,:))/ltn+lta(it-1,:);
        fdata(it,:)=sta(it,:)./lta(it,:);
    end
    fdata(1:stn+ltn+50,:)=0;
elseif mtp==3
    % calculate the classical STA/LTA with non-overlapping time window
    for it=ltn:NT-stn+1
        sta_ave=mean(tdata(it:it+stn-1,:)); % short time average
        lta_ave=mean(tdata(it-ltn+1:it,:)); % long time average
        lta_ave(lta_ave<eps)=eps; % avoid divided by 0
        fdata(it,:)=sta_ave./lta_ave;
    end
elseif mtp==4
    % no STA/LTA calculation, output the tansferred data
    fdata=tdata;
else
    error('Wrong input for mtp!');
end

end