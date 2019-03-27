function trdis=trvel2dis(trvel,dt)
% This function is used to transform the traces from particle velocity to
% particle displacement.
% INPUT-------------------------------------------------------
% trvel: input traces of particle velocity, nt*nrec;
% dt: time sampling interval;
% OUTPUT----------------------------------------------------
% trdis: output traces of particle displacement, nt*nrec.

[nt,nrec]=size(trvel); % obtain the total number of time samplings and traces

trdis=zeros(nt,nrec);

for ii=2:nt
    trdis(ii,:)=trdis(ii-1,:)+0.5*(trvel(ii-1,:)+trvel(ii,:))*dt;
end

end