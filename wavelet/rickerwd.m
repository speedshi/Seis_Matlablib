function rstf=rickerwd(freq,dt,nt,t0)
% This function is used to generate the derivatives of Ricker wavelet using
% analytical expressions of the Ricker funciton.
% freq: main frequency of Ricker wavelet (Hz);
% dt: interval of time samples (s);
% nt: number of time samples;
% t0: the input time delay (s) (default value: 0s).
% The time delay of Ricker wavelet is automatically set as 1.1/freq+t0.
% The first time sample point is set as 0s.

if nargin<4
    t0=0;
end

td=1.1/freq+t0;
factor=pi*freq;
rstf=zeros(nt,1);
for i=0:nt-1
    ft=dt*i-td;    
    rstf(i+1)=(-6*factor^2*ft+4*factor^4*ft^3)*exp(-factor^2*ft^2);
end

end