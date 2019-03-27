function [fdata1,fdata2,fdata3]=gsynwhomo_rickerw(vp,vs,den,recp,soup,mt,dt,Nt,freq)
% This function is used to generate the synthetic seismic data at a receiver point using Green's funciton.
% All parameters use SI unit, i.e. meter-kg-second.
% The wavelet is set default as analytical Ricker wavelet.
% Note: for Ricker wavelet, default time delay is 1.1/f, should keep consistent.
% All the calculation related about wavelet is performed using analytical Ricker wavelet equation.
% We have to specify the source point and the receiver point.
% Coordinate: 1-X  2-Y  3-Z
% theta: zenith angle measured form the positive Z-axis, [0, pi]
% phi: azimuth angle measured counterclockwise form the positive X-axis, [0, 2pi)
% theta & phi is used to construct direction cosines of source-receiver vector
% Input parameters: ---------------------------------------------------
% Elastic porperty of the medium:
% vp: P-wave velocity (m/s)
% vs: S-wave velocity (m/s)
% den: density (kg/m^3)
% Layout of the geometry:
% recp: position of the receivers (m) (dimension: Nre*3, Nre is the number of receivers, 3 is the X-Y-Z components in cartesion coordinate)
% soup: position of thes source point (m) (dimension: 1*3, X-Y-Z)
% Modeling parameters:
% mt: source moment tensor (dimension: 3*3), unit: N*m
% dt: time interval
% Nt: number of time samples
% freq: main frequency of the wavelet
% Output parameters: ------------------------------------------------
% (dimension: Nt*Nre)
% fdata1: displacement in X direction, unit: meter
% fdata2: displacement in Y direction, unit: meter
% fdata3: displacement in Z direction, unit: meter

Nre=size(recp,1);
piden=4*pi*den;

fdata1=zeros(Nt,Nre);
fdata2=zeros(Nt,Nre);
fdata3=zeros(Nt,Nre);

for ire=1:Nre
    rsx=recp(ire,1)-soup(1); % distance between source and receiver point in X direction
    rsy=recp(ire,2)-soup(2); % distance between source and receiver point in Y direction
    rsz=recp(ire,3)-soup(3); % distance between source and receiver point in Z direction
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
        
    tap=rd/vp; % P-wave arrival-time (s)
    tas=rd/vs; % S-wave arrival-time (s)
    
    % calculate time-dependent term
    % for near-field
    stfn=zeros(Nt,1);
    for int=1:Nt
        ingt=dt*(int-1); % time for integration
        stfn(int)=calnint(ingt,freq,tap,tas);
    end
    % for intermediate-field P-wave
    stfip=rickerw(freq,dt,Nt,tap);
    % for intermediate-field S-wave
    stfis=rickerw(freq,dt,Nt,tas);
    % for far-field P-wave
    stffp=rickerwd(freq,dt,Nt,tap);
    % for far-field S-wave
    stffs=rickerwd(freq,dt,Nt,tas);
    
    % displacement field in the near-field
    fnc1=anc1*stfn/(piden*rd^4);
    fnc2=anc2*stfn/(piden*rd^4);
    fnc3=anc3*stfn/(piden*rd^4);
    % P-wave displacement field in the intermediate-field
    fip1=aip1*stfip/(piden*vp^2*rd^2);
    fip2=aip2*stfip/(piden*vp^2*rd^2);
    fip3=aip3*stfip/(piden*vp^2*rd^2);
    % S-wave displacement field in the intermediate-field
    fis1=-ais1*stfis/(piden*vs^2*rd^2);
    fis2=-ais2*stfis/(piden*vs^2*rd^2);
    fis3=-ais3*stfis/(piden*vs^2*rd^2);
    % P-wave displacement field in the far-field
    ffp1=afp1*stffp/(piden*vp^3*rd);
    ffp2=afp2*stffp/(piden*vp^3*rd);
    ffp3=afp3*stffp/(piden*vp^3*rd);
    % S-wave displacement field in the far-field
    ffs1=-afs1*stffs/(piden*vs^3*rd);
    ffs2=-afs2*stffs/(piden*vs^3*rd);
    ffs3=-afs3*stffs/(piden*vs^3*rd);
    
    % calculate the total displacement field
    fdata1(:,ire)=fnc1+fip1+fis1+ffp1+ffs1;
    fdata2(:,ire)=fnc2+fip2+fis2+ffp2+ffs2;
    fdata3(:,ire)=fnc3+fip3+fis3+ffp3+ffs3;
end

end