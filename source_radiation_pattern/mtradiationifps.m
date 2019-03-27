function mtradiationifps(mti)
% this function is used to calculate and plot the intermediate-field S-wave and far-field P- and S-wave
% radiation pattern of moment tensor source. (S-wave is the total S-wave.)
% mti is the input moment tensor (3*3  symmetric)
% coordinate: 1-X  2-Y  3-Z
% radiation pattern based on Aki & Richards Equation 4.96.

% normalize the input moment tensor
mt=mtnorm(mti);

% create grids
nn=500;
mm=500;
theta=linspace(0,pi,nn);% zenith angle measured form the positive Z-axis, [0, pi]
phi=linspace(0,2*pi,mm);% azimuth angle measured from the positive X-axis to positive Y-axis, [0, 2pi)
[gtheta,gphi]=meshgrid(theta,phi);
[nr,nc]=size(gtheta);
% define the color map value of the intermediate-field S-wave
isr1=zeros(nr,nc);
isr2=zeros(nr,nc);
isr3=zeros(nr,nc);
ism=zeros(nr,nc);
% define the three components of the far-field S-wave
fsr1=zeros(nr,nc);
fsr2=zeros(nr,nc);
fsr3=zeros(nr,nc);
fsm=zeros(nr,nc);

% calculate the radiation pattern
% calculate the radiation pattern of far-field P-wave
fpr=(mt(1,1)*cos(gphi).^2+mt(2,2)*sin(gphi).^2+mt(1,2)*sin(2*gphi)).*sin(gtheta).^2+mt(3,3)*cos(gtheta).^2+(mt(1,3)*cos(gphi)+mt(2,3)*sin(gphi)).*sin(2*gtheta);

for jj=1:nc
    for ii=1:nr
        % direction cosine
        gmx=sin(gtheta(ii,jj))*cos(gphi(ii,jj));
        gmy=sin(gtheta(ii,jj))*sin(gphi(ii,jj));
        gmz=cos(gtheta(ii,jj));
        gama=[gmx; gmy; gmz];
        
        for iq=1:3
            for ip=1:3
                % calculate the three components of intermediate-field S-wave
                isr1(ii,jj)=isr1(ii,jj)-(6*gama(1)*gama(ip)*gama(iq)-gama(1)*deltam(ip,iq)-gama(ip)*deltam(1,iq)-2*gama(iq)*deltam(1,ip))*mt(ip,iq);
                isr2(ii,jj)=isr2(ii,jj)-(6*gama(2)*gama(ip)*gama(iq)-gama(2)*deltam(ip,iq)-gama(ip)*deltam(2,iq)-2*gama(iq)*deltam(2,ip))*mt(ip,iq);
                isr3(ii,jj)=isr3(ii,jj)-(6*gama(3)*gama(ip)*gama(iq)-gama(3)*deltam(ip,iq)-gama(ip)*deltam(3,iq)-2*gama(iq)*deltam(3,ip))*mt(ip,iq);
                ism(ii,jj)=sqrt(isr1(ii,jj)*isr1(ii,jj)+isr2(ii,jj)*isr2(ii,jj)+isr3(ii,jj)*isr3(ii,jj));
                % calculate the three components of far-field S-wave
                fsr1(ii,jj)=fsr1(ii,jj)-(gama(1)*gama(ip)-deltam(1,ip))*gama(iq)*mt(ip,iq);
                fsr2(ii,jj)=fsr2(ii,jj)-(gama(2)*gama(ip)-deltam(2,ip))*gama(iq)*mt(ip,iq);
                fsr3(ii,jj)=fsr3(ii,jj)-(gama(3)*gama(ip)-deltam(3,ip))*gama(iq)*mt(ip,iq);
                fsm(ii,jj)=sqrt(fsr1(ii,jj)*fsr1(ii,jj)+fsr2(ii,jj)*fsr2(ii,jj)+fsr3(ii,jj)*fsr3(ii,jj));
            end
        end
        
    end
end

% transfor to cartesian coordinate
% note the defined difference in matlab about polar angle 'theta'
[fpx,fpy,fpz]=sph2cart(gphi,pi/2-gtheta,abs(fpr));
[isx,isy,isz]=sph2cart(gphi,pi/2-gtheta,ism);
[fsx,fsy,fsz]=sph2cart(gphi,pi/2-gtheta,fsm);

% plot radiation pattern
% intermediate-field S-wave
figure; his=surf(isx,isy,isz,ism);
load('cmapmtrdp2','cmapmtrdp2');% 256 colormap
colormap(gca,cmapmtrdp2);
set(his,'LineStyle','none');
set(his,'FaceColor','interp');
axis equal;
box on;ax=gca;ax.BoxStyle='full';
xlabel('X');ylabel('Y');zlabel('Z');
title('intermediate-field S-wave Polarization Pattern');
mtrdpfax(isx,isy,isz,ism);% create figure which has specific axis

% far-field P-wave
figure; hp=surf(fpx,fpy,fpz,fpr);
if min(min(fpr))>=0
    % if polarization is all positive (for example an explosive source), use a specific colormap
    load('cmapmtrdpos2','cmapmtrdpos2');
    colormap(gca,cmapmtrdpos2);
elseif max(max(fpr))<=0
    % if polarization is all negtive (for example an implosion source), use a specific colormap
    load('cmapmtrdneg','cmapmtrdneg');
    colormap(gca,cmapmtrdneg);
else
    % have both positive and negtive polarization in different direction
    load('cmapmtrdp2','cmapmtrdp2');% 256 colormap
    colormap(gca,cmapmtrdp2);
end
set(hp,'LineStyle','none');
set(hp,'FaceColor','interp');
axis equal;
box on;ax=gca;ax.BoxStyle='full';
xlabel('X');ylabel('Y');zlabel('Z');
title('far-field P-wave Polarization Pattern');
mtrdpfax(fpx,fpy,fpz,fpr);% create figure which has specific axis

% far-field S-wave
figure; hfs=surf(fsx,fsy,fsz,fsm);
load('cmapmtrdp2','cmapmtrdp2');% 256 colormap
colormap(gca,cmapmtrdp2);
set(hfs,'LineStyle','none');
set(hfs,'FaceColor','interp');
axis equal;
box on;ax=gca;ax.BoxStyle='full';
xlabel('X');ylabel('Y');zlabel('Z');
title('far-field S-wave Polarization Pattern');
mtrdpfas(fsx,fsy,fsz,fsm);% create figure which has specific axis

end