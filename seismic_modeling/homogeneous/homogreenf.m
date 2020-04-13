function [ux,uy,uz]=homogreenf(vp,vs,den,recp,soup,mt,dt,Nt,t0)
% This function is used to calculate the Green's function in the
% homogeneous isotropic medium.
% All parameters use SI unit, i.e. meter-kg-second.
% Note: for Green's function the source time function is impulse funciton.
% The time serials range from: [0 , Nt-1]*dt.
% Coordinate specification: 1-X  2-Y  3-Z
% theta: zenith angle measured form the positive Z-axis, [0, pi]
% phi: azimuth angle measured from the positive X-axis to positive Y-axis, [0, 2pi)
% theta & phi is used to construct direction cosines of source-receiver vector
% INPUT--------------------------------------------------
% vp: P-wave velocity (m/s), scalar;
% vs: S-wave velocity (m/s), scalar;
% den: density (kg/m^3), scalar;
% recp: position of the receivers (m) (dimension: Nre*3, Nre is the number
% of receivers, 3 is the X-Y-Z components in cartesion coordinate);
% soup: position of thes source point (m) (dimension: 1*3, X-Y-Z);
% mt: source moment tensor (dimension: 3*3), unit: N*m;
% dt: time interval, in second;
% Nt: number of time samples;
% t0: origin time or delayed time of the source, in second.
% Output parameters:------------------------------------------
% (dimension: Nt*Nre)
% ux: displacement in X direction, unit: meter
% uy: displacement in Y direction, unit: meter
% uz: displacement in Z direction, unit: meter

Nre=size(recp,1);
piden=4*pi*den;

ux=zeros(Nt,Nre);
uy=zeros(Nt,Nre);
uz=zeros(Nt,Nre);

diamp=1/(2*dt); % the amplitude of the derivative of the impluse function

for ire=1:Nre
    rsx=recp(ire,1)-soup(1); % distance between source and receiver point in X direction
    rsy=recp(ire,2)-soup(2); % distance between source and receiver point in Y direction
    rsz=recp(ire,3)-soup(3);  % distance between source and receiver point in Z direction
    rsxy=sqrt(rsx*rsx+rsy*rsy); % distance between source and receiver point in horizontal direction (X-Y plane)
    rd=sqrt(rsx*rsx+rsy*rsy+rsz*rsz); % distance between source and receiver point
    stheta=rsxy/rd; % sin(theta)
    ctheta=rsz/rd;  % cos(theta)
    if abs(rsxy)>eps
        sphi=rsy/rsxy; % sin(phi)
        cphi=rsx/rsxy; % cos(phi)
    else
        % rsxy=0, receiver point is just above or blow the source point
        % such that sin(phi)=cos(phi)=0
        sphi=0;
        cphi=0;
    end
    % construct direction cosines
    gama(1)=stheta*cphi;
    gama(2)=stheta*sphi;
    gama(3)=ctheta;
    
    anc1=0; anc2=0; anc3=0; % near-field radiation pattern
    aip1=0; aip2=0; aip3=0; % intermediate-field P-wave radiation pattern
    ais1=0; ais2=0; ais3=0; % intermediate-field S-wave radiation pattern
    afp1=0; afp2=0; afp3=0; % far-field P-wave radiation pattern
    afs1=0; afs2=0; afs3=0; % far-field S-wave radiation pattern
    for ip=1:3
        for iq=1:3
            % calculate near-field term
            anc1=anc1+(15*gama(1)*gama(ip)*gama(iq)-3*gama(1)*deltam(ip,iq)-3*gama(ip)*deltam(1,iq)-3*gama(iq)*deltam(1,ip))*mt(ip,iq);
            anc2=anc2+(15*gama(2)*gama(ip)*gama(iq)-3*gama(2)*deltam(ip,iq)-3*gama(ip)*deltam(2,iq)-3*gama(iq)*deltam(2,ip))*mt(ip,iq);
            anc3=anc3+(15*gama(3)*gama(ip)*gama(iq)-3*gama(3)*deltam(ip,iq)-3*gama(ip)*deltam(3,iq)-3*gama(iq)*deltam(3,ip))*mt(ip,iq);
            % calculate intermediate-field P-wave term
            aip1=aip1+(6*gama(1)*gama(ip)*gama(iq)-gama(1)*deltam(ip,iq)-gama(ip)*deltam(1,iq)-gama(iq)*deltam(1,ip))*mt(ip,iq);
            aip2=aip2+(6*gama(2)*gama(ip)*gama(iq)-gama(2)*deltam(ip,iq)-gama(ip)*deltam(2,iq)-gama(iq)*deltam(2,ip))*mt(ip,iq);
            aip3=aip3+(6*gama(3)*gama(ip)*gama(iq)-gama(3)*deltam(ip,iq)-gama(ip)*deltam(3,iq)-gama(iq)*deltam(3,ip))*mt(ip,iq);
            % calculate intermediate-field S-wave term
            ais1=ais1+(6*gama(1)*gama(ip)*gama(iq)-gama(1)*deltam(ip,iq)-gama(ip)*deltam(1,iq)-2*gama(iq)*deltam(1,ip))*mt(ip,iq);
            ais2=ais2+(6*gama(2)*gama(ip)*gama(iq)-gama(2)*deltam(ip,iq)-gama(ip)*deltam(2,iq)-2*gama(iq)*deltam(2,ip))*mt(ip,iq);
            ais3=ais3+(6*gama(3)*gama(ip)*gama(iq)-gama(3)*deltam(ip,iq)-gama(ip)*deltam(3,iq)-2*gama(iq)*deltam(3,ip))*mt(ip,iq);
            % calculate far-field P-wave term
            afp1=afp1+gama(1)*gama(ip)*gama(iq)*mt(ip,iq);
            afp2=afp2+gama(2)*gama(ip)*gama(iq)*mt(ip,iq);
            afp3=afp3+gama(3)*gama(ip)*gama(iq)*mt(ip,iq);
            % calculate far-field S-wave term
            afs1=afs1+(gama(1)*gama(ip)-deltam(1,ip))*gama(iq)*mt(ip,iq);
            afs2=afs2+(gama(2)*gama(ip)-deltam(2,ip))*gama(iq)*mt(ip,iq);
            afs3=afs3+(gama(3)*gama(ip)-deltam(3,ip))*gama(iq)*mt(ip,iq);
        end
    end
    
    % calculate the time-dependent term
    tap=rd/vp; % P-wave traveltime (s)
    tas=rd/vs; % S-wave traveltime (s)
    ntp=round((t0+tap)/dt)+1; % arrival-time points for P-wave, ntp=1 means time 0.
    nts=round((t0+tas)/dt)+1; % arrival-time points for S-wave
    
    % for near-field term
    for itt=ntp:nts
        % Note: when generating Green's function, additional 'dt' is used here for 
        % convolving with source time function later (when generate synthetic data
        % using Green's function) to calculate the integration in the near-field term.
        ux(itt,ire)=ux(itt,ire)+anc1*(itt-1)*dt*dt/(piden*rd^4);
        uy(itt,ire)=uy(itt,ire)+anc2*(itt-1)*dt*dt/(piden*rd^4);
        uz(itt,ire)=uz(itt,ire)+anc3*(itt-1)*dt*dt/(piden*rd^4);
    end
    
    % for intermediate-filed terms
    % P-wave
    ux(ntp,ire)=ux(ntp,ire)+aip1/(piden*vp^2*rd^2);
    uy(ntp,ire)=uy(ntp,ire)+aip2/(piden*vp^2*rd^2);
    uz(ntp,ire)=uz(ntp,ire)+aip3/(piden*vp^2*rd^2);
    % S-wave
    ux(nts,ire)=ux(nts,ire)-ais1/(piden*vs^2*rd^2);
    uy(nts,ire)=uy(nts,ire)-ais2/(piden*vs^2*rd^2);
    uz(nts,ire)=uz(nts,ire)-ais3/(piden*vs^2*rd^2);
    
    % for far-field terms    
    % P-wave
    ux(ntp-1,ire)=ux(ntp-1,ire)+afp1*diamp/(piden*vp^3*rd);
    ux(ntp+1,ire)=ux(ntp+1,ire)+afp1*(-diamp)/(piden*vp^3*rd);
    uy(ntp-1,ire)=uy(ntp-1,ire)+afp2*diamp/(piden*vp^3*rd);
    uy(ntp+1,ire)=uy(ntp+1,ire)+afp2*(-diamp)/(piden*vp^3*rd);
    uz(ntp-1,ire)=uz(ntp-1,ire)+afp3*diamp/(piden*vp^3*rd);
    uz(ntp+1,ire)=uz(ntp+1,ire)+afp3*(-diamp)/(piden*vp^3*rd);
    % S-wave
    ux(nts-1,ire)=ux(nts-1,ire)-afs1*diamp/(piden*vs^3*rd);
    ux(nts+1,ire)=ux(nts+1,ire)-afs1*(-diamp)/(piden*vs^3*rd);
    uy(nts-1,ire)=uy(nts-1,ire)-afs2*diamp/(piden*vs^3*rd);
    uy(nts+1,ire)=uy(nts+1,ire)-afs2*(-diamp)/(piden*vs^3*rd);
    uz(nts-1,ire)=uz(nts-1,ire)-afs3*diamp/(piden*vs^3*rd);
    uz(nts+1,ire)=uz(nts+1,ire)-afs3*(-diamp)/(piden*vs^3*rd);    
end


end