function data = datatransf(data,method)
% This function is used to transfer data.
%
% INPUT:-------------------------------------------------------------------
% data: original data volume;
% method: scaler, specify the method to transfer the data volume
% if method = 0: original;
% if method = 1: absolute value;
% if method = 2: non-negtive value;
% if method = 3: square value;
% if method = 4: square root;
% if method = 5: exponential value;
% if method = 6: logarithm value;
%
% OUTPUT:------------------------------------------------------------------
% data: data volume after transfering;


if nargin < 2
    method = 0;
end


switch method
    case 0
        % original data
        fprintf('Show original data.\n');
    case 1
        % absolute
        fprintf('Show absolute value of original data.\n');
        data = abs(data);
    case 2
        % non-negtive
        fprintf('Show non-negtive value of original data.\n');
        data(data<0) = 0;
    case 3
        % square
        fprintf('Show square value of original data.\n');
        data = data.^2;
    case 4
        % square root
        fprintf('Show square root of original data.\n');
        data = sqrt(data);
    case 5
        % exponential value
        fprintf('Show exponential value of original data.\n');
        data = exp(data);
    case 6
        % logarithm value
        fprintf('Show logarithm value of original data.\n');
        data = log(data);
    otherwise
        error('Incorrect input for transformation method!');
end



end