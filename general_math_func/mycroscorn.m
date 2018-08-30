function cof=mycroscorn(A)
% This function is used to calculate the the normalized cross-correlation
% matrix of an input matrix. The cross-correlation matrix is normalized by
% the norm of the traces.
% Input:----------------------------------------------------
% A: trace matrix, each column is a trace.

sdv=sqrt(sum(A.^2)); % calcualte the norm of the input matrix along each column
cof=(A'*A)./(sdv'*sdv); % calculate the normalized cross-correlation matrix

end