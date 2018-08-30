function rstf=rickerwi(freq,dt,nt,t0)
% This function is used to calculate the integrations of the Ricker wavelet
% using numerical integration.
% freq: main frequency of Ricker wavelet (Hz);
% dt: interval of time sample (s);
% nt: number of time samples;
% t0: the input time delay (s) (default value: 0s).
% The time delay of Ricker wavelet is automatically set as 1.1/freq+t0.
% The first time sample point is set as 0s.

if nargin<4
    t0=0;
end

pif=pi*freq;
td=1.1/freq+t0; % apply time dalay, must keep consistent with 'rickerw.m'
fain=@(t)((1-2*(pif*(t-td)).^2).*exp(-(pif*(t-td)).^2)); % define the analytical equation of Ricker wavelet

rstf=zeros(nt,1);
for it=1:nt
    tn=(it-1)*dt;
    rstf(it)=integral(fain,0,tn); % apply numerical integration method
end

end