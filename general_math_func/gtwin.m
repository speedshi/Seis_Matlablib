function wind=gtwin(N,tp)
% This function is used to generate a time window which has a length of
% 2*N+1. This time window is used as weithting factor. Thus the summation
% of  these weithting factors must be 1.
% Input:
% 1. N: used to specify the length of time window (2*N+1);
% 2. tp: specify the type of the time window. 0 for Hanning window, 1 for
% 1/(r+1) type; 2 for 1/(r^2+1) type.
% Output:
% wind: time window (or we can say weighting factors).

if nargin<2
    tp=0;
end

% generate the time window
if tp==0
    wind=hanning(2*N+1);
elseif tp==1
    t=-N:1:N;
    wind=1./(abs(t)+1);
elseif tp==2
    t=-N:1:N;
    wind=1./(t.*t+1);
end

wind=wind(:);
% normalize the coefficient
tnorm=sum(wind);
wind=wind/tnorm;

end