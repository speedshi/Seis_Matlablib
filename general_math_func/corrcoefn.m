function cof=corrcoefn(A)
% This function is used to calculate the multidimensional correlation
% coefficient of the input matrix A.
% Input:----------------------------------------------------
% A: trace matrix, each column is a trace (dimension: NT*NREC).
% Output:--------------------------------------------------
% cof: a scale, is the NREC-dimensional correlation coefficient of the
% input matrix A.

N=size(A,1)-1; % used to calculate covariance N=NT-1
sdv=std(A); % calculate the standard deviation of each trace
mA=bsxfun(@minus,A,mean(A)); % remove mean value for each trace
cof=sum(prod(mA,2))/N/prod(sdv); % calculate the NREC-dimensional correlation coefficient

end