function ccv=corrcoefnv(A,n)
% This function is used to calculate all the n-dimensional correlation
% coefficients of a input matix A.
% Input:----------------------------------------------------
% A: trace matrix, each column is a trace (dimension: NT*NREC).
% n: the dimension of cross-correlation, >=2; note if n=2, use 'mycorrcoef' function is better.
% Output:--------------------------------------------------
% ccv: a vector which contains all the n-dimensional correlation
% coefficients of the input matix A.

NREC=size(A,2); % number of receivers
if NREC<n
    error('Fatal error! The input dimension n should be smaller than the number of traces.');
end

ncomb=nchoosek(1:NREC,n); % all the possible combinations of station groups, which contain n stations.
ncb=size(ncomb,1); % the number of all possible combinations
ccv=zeros(ncb,1);
for ii=1:ncb
    ccv(ii)=corrcoefn(A(:,ncomb(ii,:)));
end

end