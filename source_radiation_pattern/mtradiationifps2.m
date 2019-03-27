function mtradiationifps2(mti,azi)
% this function is used to calculate and plot the intermediate-field S-wave and far-field P- and S-wave
% radiation pattern of moment tensor source in 2-dimensional (which is an
% azimuth slice of 3D version).
% mti is the input moment tensor (3*3  symmetric)
% azi is the azimuth slice that you want to plot (in degree)
% coordinate: 1-X  2-Y  3-Z
% radiation pattern based on Aki & Richards Equation 4.96.

% normalize the input moment tensor
mt=mtnorm(mti);
azi1=azi*pi/180;
azi2=azi1+pi;

% create grids
mm=500;
nn=2*mm;
gtheta=zeros(nn,1);
gtheta(1:mm)=linspace(0,pi,mm);% zenith angle measured form the positive Z-axis, [0, pi]
gtheta(mm+1:nn)=linspace(pi,0,mm);% zenith angle measured form the positive Z-axis, [0, pi]
gphi=zeros(nn,1);% azimuth angle measured from the positive X-axis to positive Y-axis, [0, 2pi)
gphi(1:mm)=azi1;% azimuth angle measured from the positive X-axis to positive Y-axis, [0, 2pi)
gphi(mm+1:nn)=azi2;% azimuth angle measured from the positive X-axis to positive Y-axis, [0, 2pi)

% define the color map value of the intermediate-field S-wave
isr1=zeros(nn,1);
isr2=zeros(nn,1);
isr3=zeros(nn,1);
ism=zeros(nn,1);
% define the three components of the far-field S-wave
fsr1=zeros(nn,1);
fsr2=zeros(nn,1);
fsr3=zeros(nn,1);
fsm=zeros(nn,1);

% calculate the radiation pattern
% calculate the radiation pattern of far-field P-wave
fpr=(mt(1,1)*cos(gphi).^2+mt(2,2)*sin(gphi).^2+mt(1,2)*sin(2*gphi)).*sin(gtheta).^2+mt(3,3)*cos(gtheta).^2+(mt(1,3)*cos(gphi)+mt(2,3)*sin(gphi)).*sin(2*gtheta);


for ii=1:nn
    % direction cosine
    gmx=sin(gtheta(ii))*cos(gphi(ii));
    gmy=sin(gtheta(ii))*sin(gphi(ii));
    gmz=cos(gtheta(ii));
    gama=[gmx; gmy; gmz];
    
    for iq=1:3
        for ip=1:3
            % calculate the three components of intermediate-field S-wave
            isr1(ii)=isr1(ii)-(6*gama(1)*gama(ip)*gama(iq)-gama(1)*deltam(ip,iq)-gama(ip)*deltam(1,iq)-2*gama(iq)*deltam(1,ip))*mt(ip,iq);
            isr2(ii)=isr2(ii)-(6*gama(2)*gama(ip)*gama(iq)-gama(2)*deltam(ip,iq)-gama(ip)*deltam(2,iq)-2*gama(iq)*deltam(2,ip))*mt(ip,iq);
            isr3(ii)=isr3(ii)-(6*gama(3)*gama(ip)*gama(iq)-gama(3)*deltam(ip,iq)-gama(ip)*deltam(3,iq)-2*gama(iq)*deltam(3,ip))*mt(ip,iq);
            ism(ii)=sqrt(isr1(ii)*isr1(ii)+isr2(ii)*isr2(ii)+isr3(ii)*isr3(ii));
            % calculate the three components of far-field S-wave
            fsr1(ii)=fsr1(ii)-(gama(1)*gama(ip)-deltam(1,ip))*gama(iq)*mt(ip,iq);
            fsr2(ii)=fsr2(ii)-(gama(2)*gama(ip)-deltam(2,ip))*gama(iq)*mt(ip,iq);
            fsr3(ii)=fsr3(ii)-(gama(3)*gama(ip)-deltam(3,ip))*gama(iq)*mt(ip,iq);
            fsm(ii)=sqrt(fsr1(ii)*fsr1(ii)+fsr2(ii)*fsr2(ii)+fsr3(ii)*fsr3(ii));
        end
    end
    
end


% transfor to cartesian coordinate
% note the defined difference in matlab about polar angle 'theta'
%[fpx,fpy,fpz]=sph2cart(gphi,pi/2-gtheta,abs(fpr));
fpx=abs(fpr).*sin(gtheta).*cos(gphi);
fpz=abs(fpr).*cos(gtheta);
%[isx,isy,isz]=sph2cart(gphi,pi/2-gtheta,ism);
isx=abs(ism).*sin(gtheta).*cos(gphi);
isz=abs(ism).*cos(gtheta);
%[fsx,fsy,fsz]=sph2cart(gphi,pi/2-gtheta,fsm);
fsx=abs(fsm).*sin(gtheta).*cos(gphi);
fsz=abs(fsm).*cos(gtheta);

% plot radiation pattern
% intermediate-field S-wave
figure;plot(isx,isz,'k','LineWidth',2); hold on;
% far-field P-wave
plot(fpx,fpz,'r','LineWidth',2); hold on;
% far-field S-wave
plot(fsx,fsz,'b','LineWidth',2);hold off;
axis equal;
legend('intermediate S-wave','far P-wave','far S-wave')


end