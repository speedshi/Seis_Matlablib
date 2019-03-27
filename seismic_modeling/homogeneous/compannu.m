% This program is used to compare the integration results (which is required in Aki & Richards Equation 4.29 near-field term)
% between the analytic-wavelet method (calnint.m) and the numerical-wavelet method (calninum.m).
% The larger the 'r' is, the better the result is (difference is smaller).
% However using the interpolated source time function, the integration result is quite good.

clear;

vp=3500; % P-wave velocity (m/s)
vs=2000; % S-wave velocity (m/s)
dt=0.001; % time interval (s)
nt=1000; % number of time samples
freq=10; % main frequency of the analytic wavelet (Hz)
r=5; % distance between source and receiver point (m)

tlow=r/vp; % P-wave arrival-time (s)
tup=r/vs; % S-wave arrival-time (s)

stf=rickerw(freq,dt,nt); % calculate source time function

% calculate the integration using analytic wavelet
hfa=zeros(nt,1);
for ii=1:nt
    tip=(ii-1)*dt; % time where we calculate the integration, starts form 0 s.
    hfa(ii)=calnint(tip,freq,tlow,tup);
end

% calculate the integration using discretized wavelet
hfn=zeros(nt,1);
for ii=1:nt
    hfn(ii)=calninum(stf,ii,dt,tlow,tup);
end

% interpolate the source time function
N=1000; % make 'dt' become 'dt/N'
[Pstf,Pdt]=wavlintp(stf,N,dt);
% calculate the integration using discretized wavelet after interpolation
hfn2=zeros(nt,1);
for ii=1:nt
    gi=N*(ii-1)+1; % new time point for integration after interpolation
    hfn2(ii)=calninum(Pstf,gi,Pdt,tlow,tup);
end

% plot and compare
figure;plot(dt*(0:nt-1),hfa,'k','LineWidth',1.5); hold on;
plot(dt*(0:nt-1),hfn,'r','LineWidth',1.5); hold on;
plot(dt*(0:nt-1),hfn2,'b','LineWidth',1.5); hold off;
xlabel('time (s)'); ylabel('Amplitude');
legend('Analytic integration','Numerical without interpolation','Numerical with interpolation');