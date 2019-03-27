function mtradiationp(mti)
% this function is used to calculate and plot the far-field P-, S- , SV- and SH-wave
% radiation pattern of moment tensor source. (S-wave is the total S-wave.)
% mti is the input moment tensor (3*3  symmetric)
% coordinate: 1-X  2-Y  3-Z
% radiation pattern based on (C. H. Chapman 2004 Equation 4.6.18-4.6.20) or
% (Aki & Richards Equation 4.96).

% normalize the input moment tensor
mt=mtnorm(mti);

% create grids
nn=500;
mm=500;
theta=linspace(0,pi,nn);% zenith angle measured form the positive Z-axis, [0, pi]
phi=linspace(0,2*pi,mm);% azimuth angle measured from the positive X-axis to positive Y-axis, [0, 2pi)
[gtheta,gphi]=meshgrid(theta,phi);
[nr,nc]=size(gtheta);
% define the three components of the far-field S-wave
fsr1=zeros(nr,nc);
fsr2=zeros(nr,nc);
fsr3=zeros(nr,nc);
% define the color map value of the far-field S-wave
csm=zeros(nr,nc);

% calculate the radiation pattern of far-field P-wave
fpr=(mt(1,1)*cos(gphi).^2+mt(2,2)*sin(gphi).^2+mt(1,2)*sin(2*gphi)).*sin(gtheta).^2+mt(3,3)*cos(gtheta).^2+(mt(1,3)*cos(gphi)+mt(2,3)*sin(gphi)).*sin(2*gtheta);
% calculate the radiation pattern of far-field SV-wave
fsvr=0.5*(mt(1,1)*cos(gphi).^2+mt(2,2)*sin(gphi).^2-mt(3,3)+mt(1,2)*sin(2*gphi)).*sin(2*gtheta)+(mt(1,3)*cos(gphi)+mt(2,3)*sin(gphi)).*cos(2*gtheta);
% calculate the radiation pattern of far-field SH-wave
fshr=(0.5*(mt(2,2)-mt(1,1))*sin(2*gphi)+mt(1,2)*cos(2*gphi)).*sin(gtheta)+(mt(2,3)*cos(gphi)-mt(1,3)*sin(gphi)).*cos(gtheta);

% calculate the radiation pattern of far-field S-wave
for jj=1:nc
    for ii=1:nr
        % direction cosine
        gmx=sin(gtheta(ii,jj))*cos(gphi(ii,jj));
        gmy=sin(gtheta(ii,jj))*sin(gphi(ii,jj));
        gmz=cos(gtheta(ii,jj));
        gama=[gmx; gmy; gmz];
        
        for iq=1:3
            for ip=1:3               
                % calculate the three components of far-field S-wave
                if ip==1
                    fsr1(ii,jj)=fsr1(ii,jj)+(gama(iq)-gama(1)*gama(ip)*gama(iq))*mt(ip,iq);
                else
                    fsr1(ii,jj)=fsr1(ii,jj)+(-gama(1)*gama(ip)*gama(iq))*mt(ip,iq);
                end
                if ip==2
                    fsr2(ii,jj)=fsr2(ii,jj)+(gama(iq)-gama(2)*gama(ip)*gama(iq))*mt(ip,iq);
                else
                    fsr2(ii,jj)=fsr2(ii,jj)+(-gama(2)*gama(ip)*gama(iq))*mt(ip,iq);
                end
                if ip==3
                    fsr3(ii,jj)=fsr3(ii,jj)+(gama(iq)-gama(3)*gama(ip)*gama(iq))*mt(ip,iq);
                else
                    fsr3(ii,jj)=fsr3(ii,jj)+(-gama(3)*gama(ip)*gama(iq))*mt(ip,iq);
                end
                csm(ii,jj)=sqrt(fsr1(ii,jj)*fsr1(ii,jj)+fsr2(ii,jj)*fsr2(ii,jj)+fsr3(ii,jj)*fsr3(ii,jj));                
            end
        end
        
    end
end

% transfor to cartesian coordinate
% note the defined difference in matlab about polar angle 'theta'
[px,py,pz]=sph2cart(gphi,pi/2-gtheta,abs(fpr));
[sx,sy,sz]=sph2cart(gphi,pi/2-gtheta,csm);
[svx,svy,svz]=sph2cart(gphi,pi/2-gtheta,abs(fsvr));
[shx,shy,shz]=sph2cart(gphi,pi/2-gtheta,abs(fshr));

% plot radiation pattern
% P-wave
figure; hp=surf(px,py,pz,fpr);
if min(min(fpr))>=0
    % if polarization is all positive (for example an explosive source), use a specific colormap
    load('cmapmtrdpos2','cmapmtrdpos2');
    colormap(gca,cmapmtrdpos2);
else if max(max(fpr))<=0
        % if polarization is all negtive (for example an implosion source), use a specific colormap
        load('cmapmtrdneg','cmapmtrdneg');
        colormap(gca,cmapmtrdneg);
    else
        % have both positive and negtive polarization in different direction
        load('cmapmtrdp2','cmapmtrdp2');% 256 colormap
        colormap(gca,cmapmtrdp2);
    end
end
set(hp,'LineStyle','none');
set(hp,'FaceColor','interp');
axis equal; 
box on;ax=gca;ax.BoxStyle='full';
xlabel('X');ylabel('Y');zlabel('Z');
title('P-wave Polarization Pattern');
mtrdpfax(px,py,pz,fpr);% create figure which has specific axis

% S-wave
figure; has=surf(sx,sy,sz,csm);
load('cmapmtrdp2','cmapmtrdp2');% 256 colormap
colormap(gca,cmapmtrdp2);
set(has,'LineStyle','none');
set(has,'FaceColor','interp');
axis equal; 
box on;ax=gca;ax.BoxStyle='full';
xlabel('X');ylabel('Y');zlabel('Z');
title('S-wave Polarization Pattern');
mtrdpfas(sx,sy,sz,csm);% create figure which has specific axis

% SV-wave
figure; hsv=surf(svx,svy,svz,fsvr);
if min(min(fsvr))>=0
    % if polarization is all positive, use a specific colormap
    load('cmapmtrdpos2','cmapmtrdpos2');
    colormap(gca,cmapmtrdpos2);
else if max(max(fsvr))<=0
        % if polarization is all negtive, use a specific colormap
        load('cmapmtrdneg','cmapmtrdneg');
        colormap(gca,cmapmtrdneg);
    else
        % have both positive and negtive polarization in different direction
        load('cmapmtrdp2','cmapmtrdp2');% 256 colormap
        colormap(gca,cmapmtrdp2);
    end
end
set(hsv,'LineStyle','none');
set(hsv,'FaceColor','interp');
axis equal; 
box on;ax=gca;ax.BoxStyle='full';
xlabel('X');ylabel('Y');zlabel('Z');
title('SV-wave Polarization Pattern');
mtrdpfax(svx,svy,svz,fsvr);% create figure which has specific axis

% SH-wave
figure; hsh=surf(shx,shy,shz,fshr);
if min(min(fshr))>=0
    % if polarization is all positive, use a specific colormap
    load('cmapmtrdpos2','cmapmtrdpos2');
    colormap(gca,cmapmtrdpos2);
else if max(max(fshr))<=0
        % if polarization is all negtive, use a specific colormap
        load('cmapmtrdneg','cmapmtrdneg');
        colormap(gca,cmapmtrdneg);
    else
        % have both positive and negtive polarization in different direction
        load('cmapmtrdp2','cmapmtrdp2');% 256 colormap
        colormap(gca,cmapmtrdp2);
    end
end
set(hsh,'LineStyle','none');
set(hsh,'FaceColor','interp');
axis equal; 
box on;ax=gca;ax.BoxStyle='full';
xlabel('X');ylabel('Y');zlabel('Z');
title('SH-wave Polarization Pattern');
mtrdpfax(shx,shy,shz,fshr);% create figure which has specific axis
end