function f=calnint(t,freq,t1,t2)
% This function is used to calculate the integral term (in near-field term) in Aki & Richards Equation 4.29.
% The source time function used here is the analytic Ricker wavelet.
% 't' is the input time (s);
% 'freq' is the main frequency of Ricker wavelet (Hz);
% 't1' is the lower limit of integration and 't2' is the upper limit of integration, both in seconds (s);
% 'f' is the output integral value.

pif=pi*freq;
td=t-1.1/freq; % apply time dalay, must keep consistent with 'rickerw.m'
fain=@(tao)(tao.*(1-2*(pif*(td-tao)).^2).*exp(-(pif*(td-tao)).^2));
f=integral(fain,t1,t2); % apply numerical integration method

end