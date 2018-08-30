function nmt=mtnorm(mt)
% This function is used to normalize the input moment tensor such that
% the Frobenius norm of moment tensor is 1, which means sum(m_ij^2)=1.

[nr,nc]=size(mt);
tempc=0;
for j=1:nc
    for i=1:nr
        tempc=tempc+mt(i,j)*mt(i,j);
    end
end
nmt=mt/sqrt(tempc);
end