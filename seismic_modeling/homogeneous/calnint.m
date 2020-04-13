function f=calnint(t,freq,t1,t2,t0)
% This function is used to calculate the integral term (in near-field term) in Aki & Richards Equation 4.29.
% The source time function used here is the analytic Ricker wavelet.
% Note: for Ricker wavelet, the time delay is (1.1/f+t0), should keep
% consistent with "rickerw.m".
% INPUT:
% t: the input time (s);
% freq: the main frequency of Ricker wavelet (Hz);
% t1: the lower limit of integration and 't2' is the upper limit of integration, both in seconds (s);
% t0: origin time or delayed time of the wavelet, in second.
% OUTPUT:
% f: the output integral value.

if nargin<5
   t0=0; 
end

pif=pi*freq;
td=t-1.1/freq-t0; % apply time dalay, must keep consistent with 'rickerw.m'
fain=@(tao)(tao.*(1-2*(pif*(td-tao)).^2).*exp(-(pif*(td-tao)).^2));
f=integral(fain,t1,t2); % apply numerical integration method

end