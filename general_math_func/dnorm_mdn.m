function [data_norm, trend]=dnorm_mdn(data,n)
% This function is used to normalize the data by minus the data's
% median values within a sliding window.
%
% INPUT--------------------------------------------------------------------
% data: original data to be normalized, a vector: nt*1;
% n: the number of data points in a sliding window;
%
% OUTPUT-------------------------------------------------------------------
% data_norm: data after removing median values, nt*1;
% trend: the moving median values (data trend) which are removed from the data.

% set default parameters
if nargin<2
    n=10;
end


nt=length(data); % total number of data points (time samples)


% Initialize
trend=zeros(size(data));


if n>1
    
    if n>nt
        n=nt;
    end
    
    nphf=round(n/2); % half window length, at least one point
    
    for ii=1:nt
        
        if ii<=nphf
            trend(ii)=median(data(1:2*nphf));
        elseif ii>nt-nphf
            trend(ii)=median(data((nt-2*nphf+1):nt));
        else
            trend(ii)=median(data((ii-nphf):(ii+nphf)));
        end
        
    end
    data_norm=data-trend; % remove moving median values
    
else
    data_norm=data;
end


end