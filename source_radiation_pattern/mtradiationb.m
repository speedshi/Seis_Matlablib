function mtradiationb(mti)
% this function is used to calculate and plot the far-field P-wave
% radiation pattern of moment tensor source using beach ball.
% mti is the input moment tensor (3*3  symmetric)
% coordinate: 1-X  2-Y  3-Z
% radiation pattern based on (C. H. Chapman 2004 Equation 4.6.18) or
% (Aki & Richards Equation 4.29).
% plane beach ball projection using stereographic projection (Aki & Richards figure 4.16)

% normalize the input moment tensor
mt=mtnorm(mti);

% create grids
nn=800;
mm=800;
theta=linspace(0,pi,nn);% zenith angle measured form the positive Z-axis, [0, pi]
phi=linspace(0,2*pi,mm);% azimuth angle measured from the positive X-axis to positive Y-axis, [0, 2pi)
[gtheta,gphi]=meshgrid(theta,phi);

% set a unit sphere
[nr,nc]=size(gtheta);
runit=ones(nr,nc);

% calculate the radiation pattern of far-field P-wave
fpr=(mt(1,1)*cos(gphi).^2+mt(2,2)*sin(gphi).^2+mt(1,2)*sin(2*gphi)).*sin(gtheta).^2+mt(3,3)*cos(gtheta).^2+(mt(1,3)*cos(gphi)+mt(2,3)*sin(gphi)).*sin(2*gtheta);

% transfor to cartesian coordinate
% note the defined difference in matlab about polar angle 'theta'
[px,py,pz]=sph2cart(gphi,pi/2-gtheta,runit);
figure;surf(px,py,pz,fpr);
colormap([1,1,1;0,0,0]);
axis equal;shading interp;
axis off;hold on;

alim=1.8; % set the axis limit value
% plot the three axis
% with negtive axis
%quiver3(-alim,0,0,2*alim,0,0,'Color','k','LineWidth',2,'MaxHeadSize',0.2,'AutoScale','off'); hold on;
%quiver3(0,-alim,0,0,2*alim,0,'Color','k','LineWidth',2,'MaxHeadSize',0.2,'AutoScale','off'); hold on;
%quiver3(0,0,-alim,0,0,2*alim,'Color','k','LineWidth',2,'MaxHeadSize',0.2,'AutoScale','off'); hold on;
% without negtive axis
quiver3(0,0,0,alim,0,0,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on;
quiver3(0,0,0,0,alim,0,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on;
quiver3(0,0,0,0,0,alim,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on;

% add text label for the three axis
axtxv=1.1*alim; % control the position of the text label, slightly bigger than the axis limit
text(axtxv,0,0,'X','HorizontalAlignment','center');
text(0,axtxv,0,'Y','HorizontalAlignment','center');
text(0,0,axtxv,'Z','HorizontalAlignment','center');
hold off;




% now projection to the plane beach ball
% using stereographic projection
% projection using the lower hemisphere
% create grids, as we will plot each scatter point, the grid should be denser
bn=1000;
bm=1000;
theta=linspace(pi/2,pi,bn); % zenith angle measured form the positive Z-axis, for lower hemisphere: [pi/2, pi]
phi=linspace(0,2*pi,bm); % azimuth angle measured from the positive X-axis to positive Y-axis, [0, 2pi)
[gtheta,gphi]=meshgrid(theta,phi);
upr=ones(size(phi));

% calculate the radiation pattern of far-field P-wave
fpr=(mt(1,1)*cos(gphi).^2+mt(2,2)*sin(gphi).^2+mt(1,2)*sin(2*gphi)).*sin(gtheta).^2+mt(3,3)*cos(gtheta).^2+(mt(1,3)*cos(gphi)+mt(2,3)*sin(gphi)).*sin(2*gtheta);

bbr=tan((pi-gtheta)/2); % note zenith angle need to be transferred to take-off angle
[bbx,bby]=pol2cart(gphi,bbr); % transfer to cartesian coordinate

% plot the beach ball in a plane
% XY profile, view from positive Z axis
flg=fpr>=-eps; % find the non-negative position
figure;plot(bbx(flg),bby(flg),'.k');hold on; % plot positive value -- explosive
flg=~flg; % find the negtive position
plot(bbx(flg),bby(flg),'.w');hold on; % plot negative value -- implosive
% plot the boundary of beach ball
plh=polar(phi,upr,'k');hold on;
plh.LineWidth=1.6;
% plot the axis at the corner
baxlm=0.5; % set the axis length
baxx=-1.25; baxy=-1.25; % set the original point for the new added axis
quiver(baxx,baxy,baxlm,0,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on;
quiver(baxx,baxy,0,baxlm,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on;
% add text label for the three axis
btxv=1.15*baxlm; % control the position of the text label, slightly bigger than the axis limit
text(baxx+btxv,baxy,'X','HorizontalAlignment','center');
text(baxx,baxy+btxv,'Y','HorizontalAlignment','center');
axis tight; axis equal;
axis off;

end