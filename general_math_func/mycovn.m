function cof=mycovn(A)
% This function is used to calculate the the normalized covariance matrix
% of an input matrix. The covariance matrix is normalized by the norm of
% the traces.
% Input:----------------------------------------------------
% A: trace matrix, each column is a trace.

sdv=sqrt(sum(A.^2)); % calcualte the norm of the input matrix along each column
cof=cov(A)./(sdv'*sdv); % calculate the normalized covariance matrix

end