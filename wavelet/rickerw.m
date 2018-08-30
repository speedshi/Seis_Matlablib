function rstf=rickerw(freq,dt,nt,t0)
% This function is used to generate the Ricker wavelet.
% The time delay of Ricker wavelet is automatically set as 1.1/freq+t0.
% The first time sample point is set as 0s.
% INPUT-------------------------------------------------------
% freq: main frequency of Ricker wavelet (Hz);
% dt: interval of time samples (s);
% nt: number of time samples;
% t0: the input time delay (s) (default value: 0s).
% OUTPUT---------------------------------------------------
% rstf: Ricker wavelet.


if nargin<4
    t0=0;
end

td=1.1/freq+t0;
rstf=zeros(nt,1);

for i=0:nt-1
    ft=dt*i-td;
    factor=(pi*freq*ft)*(pi*freq*ft);
    rstf(i+1)=(1-2*factor)*exp(-factor);
end

end