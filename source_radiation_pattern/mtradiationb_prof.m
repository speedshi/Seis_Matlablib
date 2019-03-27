function mtradiationb_prof(mti)
% This function is used to calculate and plot (profiles) the far-field P-wave
% radiation pattern of moment tensor source using beach ball seen from a
% perticular directions (X or Y or Z axis).
% Coordinate: 1-X  2-Y  3-Z.
% Radiation pattern based on (C. H. Chapman 2004 Equation 4.6.18) or
% (Aki & Richards Equation 4.29).
% Using Red-blue colormap. red--black, blue--white.
% Input:----------------------------------------------------
% mti: the input moment tensor source (3*3  symmetric).
% Output:--------------------------------------------------
% three profiles of the far-field P-wave radiation beach ball seen from
% different directions.


% normalize the input moment tensor
mt=mtnorm(mti);

% create grids
nn=1000;
mm=1000;
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
colormap([0,0,1;1,0,0]);%red-blue colormap; red--black, blue--white
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




% plot the YZ profile, view from negative X axis
figure;surf(px,py,pz,fpr); hold on;
bxx=zeros(nn,1);
byy=sin(linspace(0,2*pi,nn));
bzz=cos(linspace(0,2*pi,nn));
% plot the boundary of beach ball
plh=plot3(bxx,byy,bzz,'k');hold on;
plh.LineWidth=1.8;
colormap([0,0,1;1,0,0]);%red-blue colormap; red--black, blue--white
% plot the axis at the corner
baxlm=0.5; % set the axis length
baxx=-1.65; baxy=-1.65; baxz=-1.65;% set the original point for the new added axis
quiver3(baxx,baxy,baxz,0,baxlm,0,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on; % Y axis
quiver3(baxx,baxy,baxz,0,0,baxlm,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on; % Z axis
% add text label for the three axis
axtxv=1.1*baxlm; % control the position of the text label, slightly bigger than the axis limit
text(baxx,baxy+axtxv,baxz,'Y','HorizontalAlignment','center');
text(baxx,baxy,baxz+axtxv,'Z','HorizontalAlignment','center');
axis equal;shading interp;
axis off; hold off;
view(90,180); % YZ profile, view from negative X axis


% plot the XZ profile, view from positive Y axis
figure;surf(px,py,pz,fpr); hold on;
bxx=sin(linspace(0,2*pi,nn));
byy=zeros(nn,1);
bzz=cos(linspace(0,2*pi,nn));
% plot the boundary of beach ball
plh=plot3(bxx,byy,bzz,'k');hold on;
plh.LineWidth=1.8;
colormap([0,0,1;1,0,0]);%red-blue colormap; red--black, blue--white
% plot the axis at the corner
quiver3(baxx,baxy,baxz,baxlm,0,0,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on; % X axis
%quiver3(baxx,baxy,baxz,0,baxlm,0,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on; % Y axis
quiver3(baxx,baxy,baxz,0,0,baxlm,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on; % Z axis
% add text label for the three axis
text(baxx+axtxv,baxy,baxz,'X','HorizontalAlignment','center');
%text(baxx,baxy+axtxv,baxz,'Y','HorizontalAlignment','center');
text(baxx,baxy,baxz+axtxv,'Z','HorizontalAlignment','center');
axis equal;shading interp;
axis off; hold off;
view(0,180); % XZ profile, view from positive Y axis


% plot the XY profile, view from negative Z axis
figure;surf(px,py,pz,fpr); hold on;
bxx=sin(linspace(0,2*pi,nn));
byy=cos(linspace(0,2*pi,nn));
bzz=zeros(nn,1);
% plot the boundary of beach ball
plh=plot3(bxx,byy,bzz,'k');hold on;
plh.LineWidth=1.8;
colormap([0,0,1;1,0,0]);%red-blue colormap; red--black, blue--white
% plot the axis at the corner
quiver3(baxx,baxy,baxz,baxlm,0,0,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on; % X axis
quiver3(baxx,baxy,baxz,0,baxlm,0,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on; % Y axis
%quiver3(baxx,baxy,baxz,0,0,baxlm,'Color','k','LineWidth',2,'MaxHeadSize',0.4,'AutoScale','off'); hold on; % Z axis
% add text label for the three axis
text(baxx+axtxv,baxy,baxz,'X','HorizontalAlignment','center');
text(baxx,baxy+axtxv,baxz,'Y','HorizontalAlignment','center');
%text(baxx,baxy,baxz+axtxv,'Z','HorizontalAlignment','center');
axis equal;shading interp;
axis off; hold off;
view(90,-90); % XY profile, view from negative Z axis


end