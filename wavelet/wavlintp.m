function [Pstf,Pdt]=wavlintp(stf,N,dt,cs)
% This function is used to interpolate the input source time function.
% The input source time function is discretized and we don't know the analytic expression of it.
% input parameters: ----------------------------------------------
% stf: source time function (discretized, a vector: 1-nt);
% N: time interval will be 1/N of the original one after interpolatioin
% dt: time sample interval (s);
% Output parameters: -------------------------------------------
% 'Pstf': the interpolated source time function;
% 'Pdt': the new time interval after interpolation (s).

if nargin<4
    cs=0;
end

nt=max(size(stf)); % number of time samples for the input wavelet

Pdt=dt/N; % new time interval after interpolation
Pnt=N*(nt-1)+1; % new number of time samples after interpolation

tor=dt*(0:nt-1);
tpor=Pdt*(0:Pnt-1);
Pstf=interp1(tor,stf,tpor,'spline');

if cs~=0
    % plot the interpolated source time function
    figure;plot(tor,stf,'k','LineWidth',1.5); hold on;
    plot(tpor,Pstf,'r','LineWidth',1.5); hold off;
    xlabel('time (s)'); ylabel('Amplitude');
    legend('Original wavelet','Interpolated wavelet');
end

end
