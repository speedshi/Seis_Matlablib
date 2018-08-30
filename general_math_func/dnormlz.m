function ndata=dnormlz(data,n1,n2)
% This function is used to linear normalize the input data to the specified range.
% The value of output data after normalization is between [n1 n2].
% Input:----------------------------------------------------
% data: the original input data
% n1: lower limit of the output data
% n2: upper limit of the output data
% Output:--------------------------------------------------
% ndata: output data after normalization.

dmax=max(data(:));
dmin=min(data(:));
ndata=(n2-n1)/(dmax-dmin)*data+(n1*dmax-n2*dmin)/(dmax-dmin);
end