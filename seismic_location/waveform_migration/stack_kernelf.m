function tsdata=stack_kernelf(fdata,stp,N,wtp)
% This function is used to process the input data and transfer it into
% another form. The time window is applied here. The output is the
% summation of the transfered data within the selected time window.
% Usage:
% input:----------------------------------------------------
% 1 fdata: the input original data, should be a 2-dimensional data; and the
% transformation is processed along 1st dimension (e.g. the time samples is
% stored along the 1st dimension).
% 2 stp: the number used to specify the transformation type. 0- no
% transformation; 1- absolute value; 2- envelope; 3- non-negtive value.
% 3 N: related to the length of the time window (2*N+1).
% 4 wtp: used to specify the type of the time window.
% output:--------------------------------------------------
% tdata: the data after transformation.


% set default value
if nargin<2
    stp=0;
    N=0;
    wtp=0;
elseif nargin<3
    N=0;
    wtp=0;
elseif nargin<4
    wtp=0;
end

%tdata=zeros(size(fdata));
if stp==0
    % original data
    tdata=fdata;
elseif stp==1
    % absolute value
    tdata=abs(fdata);
elseif stp==2
    % envelope
    tdata=abs(hilbert(fdata));
elseif stp==3
    % non-negtive
    tdata=fdata;
    tdata(tdata<0)=0;
else
    error('The input stp is not recognised!');
end

% Generate time window, used for weighted summation in the time window
% If N=0, then only one point will be used in the summation. Thus no time
% window will be applied in this situation.
twind=gtwin(N,wtp);

% do summation in the selected time window
% e.g. weighted summation
[Nt,nre]=size(fdata); % number of time samples and geophones
tsdata=zeros(size(fdata)); % initialize the final output data

if N>0
    % upper part
    for ii=1:N
        unt=min([ii-1  N-ii]); % new half time window length
        utwd=gtwin(unt,wtp); % new time window
        for ir=1:nre
            tsdata(ii,ir)=sum(tdata(ii-unt:ii+unt,ir).*utwd);
        end
    end
    
    % down part
    for ii=Nt-N+1:Nt
        unt=min([ii-Nt+N-1  Nt-ii]); % new half time window length
        utwd=gtwin(unt,wtp); % new time window
        for ir=1:nre
            tsdata(ii,ir)=sum(tdata(ii-unt:ii+unt,ir).*utwd);
        end
    end
end

% inner part
for ii=N+1:Nt-N
    for ir=1:nre
        tsdata(ii,ir)=sum(tdata(ii-N:ii+N,ir).*twind);
    end
end

end