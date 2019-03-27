function mtradiationv(mti)
% This function is used to calculate and plot the far-field P- and S-wave
% radiation pattern of moment tensor source in Vector model.
% S-wave is the total S-wave and will not be decomposed into SV- and SH-wave.
% Coordinate: 1-X  2-Y  3-Z
% Radiation pattern based on Aki & Richards Equation 4.29.
% INPUT--------------------------------------------------
% mti: moment tensor, 3*3  symmetric matrix.
% OUTPUT-----------------------------------------------
% vector figures showing the far-field P- and S-wave radiation pattern.
% Author: Peidong Shi (c) 2018-11-22.


% normalize the input moment tensor
mt=mtnorm(mti);

% create grids
nn=32;
mm=36;
theta=linspace(0,pi,nn);% zenith angle measured form the positive Z-axis, [0, pi]
phi=linspace(0,2*pi,mm);% azimuth angle measured from the positive X-axis to positive Y-axis, [0, 2pi)
[gtheta,gphi]=meshgrid(theta,phi);

[nr,nc]=size(gtheta);
runit=ones(nr,nc);
% define the three components of the far-field P-wave
fpr1=zeros(nr,nc);
fpr2=zeros(nr,nc);
fpr3=zeros(nr,nc);
% define the color map value of the far-field P-wave
cpm=zeros(nr,nc);
% define the three components of the far-field S-wave
fsr1=zeros(nr,nc);
fsr2=zeros(nr,nc);
fsr3=zeros(nr,nc);
% define the color map value of the far-field S-wave
csm=zeros(nr,nc);

for jj=1:nc
    for ii=1:nr
        % direction cosine
        gmx=sin(gtheta(ii,jj))*cos(gphi(ii,jj));
        gmy=sin(gtheta(ii,jj))*sin(gphi(ii,jj));
        gmz=cos(gtheta(ii,jj));
        gama=[gmx; gmy; gmz];
        
        for iq=1:3
            for ip=1:3
                % calculate the three components of far-field P-wave
                fpr1(ii,jj)=fpr1(ii,jj)+gama(1)*gama(ip)*gama(iq)*mt(ip,iq);
                fpr2(ii,jj)=fpr2(ii,jj)+gama(2)*gama(ip)*gama(iq)*mt(ip,iq);
                fpr3(ii,jj)=fpr3(ii,jj)+gama(3)*gama(ip)*gama(iq)*mt(ip,iq);                
                cpm(ii,jj)=(fpr1(ii,jj)+fpr2(ii,jj)+fpr3(ii,jj))/(gama(1)+gama(2)+gama(3));
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

% generate sphere in cartesian coordinate
% note the defined difference in matlab about polar angle 'theta'
[rnx,rny,rnz]=sph2cart(gphi,pi/2-gtheta,runit);
% load a saved colormap
load('cmapmtv4','cmapmtv');

% plot the radiation pattern
% plot Vector radiation pattern
figure; % for P-wave
% setup a specific colormap
colormap(cmapmtv);
hp=quiver3c(rnx,rny,rnz,fpr1,fpr2,fpr3,cpm);
hp.AutoScaleFactor=2.6;
set(gca, 'color', 'k');% set the background color to be black
axis equal;
xlabel('X');ylabel('Y');zlabel('Z');
%title('P-wave Polarization Pattern');

figure; % for S-wave
% setup a specific colormap
colormap(cmapmtv);
hs=quiver3c(rnx,rny,rnz,fsr1,fsr2,fsr3,csm);
hs.AutoScaleFactor=1.6;
set(gca, 'color', 'k');% set the background color to be black
axis equal;
xlabel('X');ylabel('Y');zlabel('Z');
xlim([-1 1]); set(gca,'xtick',-1:0.5:1);% set the tick of X-axis
ylim([-1 1]); set(gca,'ytick',-1:0.5:1);% set the tick of Y-axis
zlim([-1 1]); set(gca,'ztick',-1:0.5:1);% set the tick of Z-axis
%title('S-wave Polarization Pattern');

end