function f=calninum(stf,tnp,dt,t1,t2)
% This function is used to calculate the integral term (in near-field term) in Aki & Richards Equation 4.29.
% The input source time function is discretized and we don't know the analytic expression of it.
% So here the numerical integration method based on the discretized wavelet is used.
% INPUT--------------------------------------------------
% stf: source time function, a vector, nt*1;
% tnp: on which time point the integration is calculated, an integer;
% dt: time sample interval (s);
% t1: the lower limit of integration, in seconds (s);
% t2: the upper limit of integration, in seconds (s).
% OUTPUT-----------------------------------------------
% f: integral value, scalar.

nt=max(size(stf));

nlow=round(t1/dt); % time point for P-wave arrival-time t1
nup=round(t2/dt); % time point for S-wave arrival-time t2

f=0;
for i=nlow:nup
    inx=tnp-i;
    if inx>=1 && inx<=nt       
        f=f+i*dt*stf(inx)*dt;
    end
end

end