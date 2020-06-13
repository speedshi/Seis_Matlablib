function data=ireplace_zeros(data)
% This function is used to replace the zero values in the data using random
% eps values.
%
% INPUT--------------------------------------------------------------------
%      data: the input data, can be a vector, matrix or multidimensional
%      array.
%


index = (data == 0); % find the index for zero values
N = sum(index,'all'); % the number of zeros in the data

data(index) = randn(N,1)*eps;


end