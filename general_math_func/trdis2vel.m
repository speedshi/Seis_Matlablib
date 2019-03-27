function trvel=trdis2vel(trdis,dt)
% This function is used to transform the traces from particle displacement to
% particle velocity.
% INPUT-------------------------------------------------------
% trdis: input traces of particle displacement, nt*nrec;
% dt: time sampling interval;
% OUTPUT----------------------------------------------------
% trvel: output traces of particle velocity, nt*nrec.

[nt,nrec]=size(trdis); % obtain the total number of time samplings and traces

trvel=zeros(nt,nrec);

for ii=2:nt-1
    trvel(ii,:)=0.5*(trdis(ii+1,:)-trdis(ii-1,:))/dt;
end

end