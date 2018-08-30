function cof=mycorrcoef(A)
% This function is used to calculate the correlation coefficient matrix of
% an input matrix. Only the upper triangular part of the output matrix 'cof'
% matters.
% Input:----------------------------------------------------
% A: trace matrix, each column is a trace.

sdv=std(A); % calcualte the standard deviation of the input matrix
cof=cov(A)./(sdv'*sdv); % calculate the correlation coefficient matrix

end