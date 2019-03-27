function fd=waveldely(stf,tdl,dt)
% This function is used to apply a time delay on an input source time funciton.
% INPUT--------------------------------------------------
% stf: input source time function, starting from time 0;
% tdl: time delay (s);
% dt: time sample interval (s).
% OUTPUT-----------------------------------------------
% fd: output source time function after time delay.

tdi=round(tdl/dt);
nt=max(size(stf));
fd=zeros(nt,1);
fd(tdi+1:nt)=stf(1:nt-tdi);

end