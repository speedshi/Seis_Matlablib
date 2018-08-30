function disp_3dslice(xr,yr,zr,xc,yc,zc,px,py,pz)
% This function is used to display the 3D slices in one figure.
% Input:----------------------------------------------------
% xr: X-axis range (km), 2*1 vector;
% yr: Y-axis range (km), 2*1 vector;
% zr: Z-axis range (km), 2*1 vector;
% xc: the X coordinate of YZ profile (km), scalar;
% yc: the Y coordinate of XZ profile (km), scalar;
% zc: the X coordinate of XY profile (km), scalar;
% px: slice in the X direction (YZ profile), ny*nz matrix, first dimension: Y; second dimension: Z;
% py: slice in the Y direction (XZ profile), nx*nz matrix, first dimension: X; second dimension: Z;
% pz: slice in the Z direction (XY profile), nx*ny matrix, first dimension: X; second dimension: Y.

[ny,nz]=size(px); % number of data points in the Z and Y directions
nx=size(py,1); % number of data points in the X direction

xx=linspace(xr(1),xr(2),nx); % X coordinates list
yy=linspace(yr(1),yr(2),ny); % Y coordinates list
zz=linspace(zr(1),zr(2),nz); % Z coordinates list

% plot the 3 slices in one figure
figure; 
[mmz,~]=meshgrid(zz,yy);
shx=surf(xc*ones(size(zz)),yy,mmz,px); hold on;
shx.LineStyle='none';

[~,mmz]=meshgrid(xx,zz);
shy=surf(xx,yc*ones(size(zz)),mmz,py'); hold on;
shy.LineStyle='none';

shz=surf(xx,yy,zc*ones(length(yy),length(xx)),pz'); hold on;
shz.LineStyle='none';

xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
% title('Envelope');
axis tight equal;
box on; grid on; %grid minor;
ax=gca; 
ax.GridLineStyle=':';
%ax.BoxStyle = 'full';
colormap('jet(512)');colorbar;

end