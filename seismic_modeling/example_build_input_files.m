% This function is used to generate the receiver and velocity files for modeling.
% The unit for distance is meter, for velocity is m/s, for density is kg/m^3, SI unit.

% Set some parameters------------------------------------------------------
xs=2500; % X coordinate of the source
ys=2500; % Y coordinate of the source
zs=1050; % Z coordinate of the source
% below is the parameters used in the 3D modelling
xm0=0; % reference origin point in the X direction, meter
ym0=0; % reference origin point in the Y direction, meter
zm0=0; % reference origin point in the Z direction, meter
lx=5000; % length of the model in the X direction, meter
ly=5000; % length of the model in the Y direction, meter
lz=3000; % length of the model in the Z direction, meter
dx=5; % grid interval in the X direction, meter
dy=5; % grid interval in the Y direction, meter
dz=5; % grid interval in the Z direction, meter
%-------------------------------------------------------------------------------------------

% Generate the receiver file for homogeneous and layered models------------
fid1=fopen('receiver1.dat','wt'); % file name of the receiver file
nn=5*8+3*51+8*2; % total number of receivers
fprintf(fid1,'%d\n',nn); % write the total number of stations into the file

% For vertical array 1:
nr=8; % number of receivers
xr=[4000 4000]; % X range of the array
yr=[1000 1000]; % Y range of the array
zr=[650   1000]; % Z range of the array
xx=linspace(xr(1),xr(2),nr);
yy=linspace(yr(1),yr(2),nr);
zz=linspace(zr(1),zr(2),nr);
% write on the receiver file
id=0;
for ii=1:nr
    id=id+1; % index of the station
    fprintf(fid1,'s%03d  %9.3f  %9.3f  %9.3f\n',id,xx(ii),yy(ii),zz(ii));
end

% For vertical array 2:
nr=8; % number of receivers
xr=[3500 3500]; % X range of the array
yr=[3500 3500]; % Y range of the array
zr=[650   1000]; % Z range of the array
xx=linspace(xr(1),xr(2),nr);
yy=linspace(yr(1),yr(2),nr);
zz=linspace(zr(1),zr(2),nr);
% write on the receiver file
for ii=1:nr
    id=id+1; % index of the station
    fprintf(fid1,'s%03d  %9.3f  %9.3f  %9.3f\n',id,xx(ii),yy(ii),zz(ii));
end

% For vertical array 3:
nr=8; % number of receivers
xr=[2500 2500]; % X range of the array
yr=[4000 4000]; % Y range of the array
zr=[650   1000]; % Z range of the array
xx=linspace(xr(1),xr(2),nr);
yy=linspace(yr(1),yr(2),nr);
zz=linspace(zr(1),zr(2),nr);
% write on the receiver file
for ii=1:nr
    id=id+1; % index of the station
    fprintf(fid1,'s%03d  %9.3f  %9.3f  %9.3f\n',id,xx(ii),yy(ii),zz(ii));
end

% For vertical array 4:
nr=8; % number of receivers
xr=[1500 1500]; % X range of the array
yr=[2500 2500]; % Y range of the array
zr=[650   1000]; % Z range of the array
xx=linspace(xr(1),xr(2),nr);
yy=linspace(yr(1),yr(2),nr);
zz=linspace(zr(1),zr(2),nr);
% write on the receiver file
for ii=1:nr
    id=id+1; % index of the station
    fprintf(fid1,'s%03d  %9.3f  %9.3f  %9.3f\n',id,xx(ii),yy(ii),zz(ii));
end

% For vertical array 5:
nr=8; % number of receivers
xr=[2000 2000]; % X range of the array
yr=[1000 1000]; % Y range of the array
zr=[650   1000]; % Z range of the array
xx=linspace(xr(1),xr(2),nr);
yy=linspace(yr(1),yr(2),nr);
zz=linspace(zr(1),zr(2),nr);
% write on the receiver file
for ii=1:nr
    id=id+1; % index of the station
    fprintf(fid1,'s%03d  %9.3f  %9.3f  %9.3f\n',id,xx(ii),yy(ii),zz(ii));
end

% For horizontal array 1:
nr=51; % number of receivers
xr=[3250 3250]; % X range of the array
yr=[1000 3500]; % Y range of the array
zr=[1040 1040]; % Z range of the array
xx=linspace(xr(1),xr(2),nr);
yy=linspace(yr(1),yr(2),nr);
zz=linspace(zr(1),zr(2),nr);
% write on the receiver file
for ii=1:nr
    id=id+1; % index of the station
    fprintf(fid1,'s%03d  %9.3f  %9.3f  %9.3f\n',id,xx(ii),yy(ii),zz(ii));
end

% For horizontal array 2:
nr=51; % number of receivers
xr=[2500 2500]; % X range of the array
yr=[1000 3500]; % Y range of the array
zr=[1040 1040]; % Z range of the array
xx=linspace(xr(1),xr(2),nr);
yy=linspace(yr(1),yr(2),nr);
zz=linspace(zr(1),zr(2),nr);
% write on the receiver file
for ii=1:nr
    id=id+1; % index of the station
    fprintf(fid1,'s%03d  %9.3f  %9.3f  %9.3f\n',id,xx(ii),yy(ii),zz(ii));
end

% For horizontal array 3:
nr=51; % number of receivers
xr=[1750 1750]; % X range of the array
yr=[1000 3500]; % Y range of the array
zr=[1040 1040]; % Z range of the array
xx=linspace(xr(1),xr(2),nr);
yy=linspace(yr(1),yr(2),nr);
zz=linspace(zr(1),zr(2),nr);
% write on the receiver file
for ii=1:nr
    id=id+1; % index of the station
    fprintf(fid1,'s%03d  %9.3f  %9.3f  %9.3f\n',id,xx(ii),yy(ii),zz(ii));
end

% For spherical array:
azi=0:45:360-45; % generate azimuth angles, measured from the positive X-axis to positive Y-axis (positive X-axis is 0, positive Y-axis is pi/2)
zen=30:30:60; % generate zenith angles, measured form the up-direction (vertical up: 0; horizontal: 90; verical down: 180), upper hemisphere
naz=length(azi); % number of azimuth angles
nze=length(zen); % number of zenith angles
nr=naz*nze; % number of receivers
dist=400; % distance between the recreiver and the source
for ia=1:naz
    for iz=1:nze
        id=id+1; % index of the station
        dxy=dist*sind(zen(iz)); % horizontal distance
        xx=xs+dxy*cosd(azi(ia)); % X coordinate of the receiver
        yy=ys+dxy*sind(azi(ia)); % Y coordinate of the receiver
        zz=zs-dist*cosd(zen(iz)); % Z coordinate of the reveiver
        % place the spherical arrays to the nearest grid points (grid point has a interval of 5 m)
        xx=round(xx/5)*5;
        yy=round(yy/5)*5;
        zz=round(zz/5)*5;
        fprintf(fid1,'s%03d  %9.3f  %9.3f  %9.3f\n',id,xx,yy,zz); % write to the file
    end
end

fclose(fid1);
%-------------------------------------------------------------------------------------------------------------


% Generate the receiver file for 3D modeling. Use grid points rather than
% the absolute positions---------------------------------------------------
fid1=fopen('receiver.dat','wt'); % file name of the receiver file
nn=dlmread('receiver1.dat','',[0 0 0 0]); % read in the total number of the receivers
fprintf(fid1,'%d\n',nn); % write the total number of stations into the file
recp=dlmread('receiver1.dat','',1,0); % read in the absolute positions of the receivers
xgd=(recp(:,1)-xm0)/dx+1; % receiver grid point in the X direction
ygd=(recp(:,2)-ym0)/dy+1; % receiver grid point in the Y direction
zgd=(recp(:,3)-zm0)/dz+1; % receiver grid point in the Z direction
for ii=1:nn
    fprintf(fid1,'%5d  %5d  %5d\n',xgd(ii),ygd(ii),zgd(ii)); % write the receiver grid points to the file
end
fclose(fid1);
%---------------------------------------------------------------------------------------------


% Generate velocity file for homogeneous model-----------------------------
vp=4000; % P-wave velocity
vs=2200; % S-wave velocity
den=2000; % density
sdep=0; % elevation of the free surface
fid=fopen('homo.dat','wt'); % file name of the velocity model
fprintf(fid,'%9.3f\n',sdep);
fprintf(fid,'%9.3f  %9.3f  %9.3f %9.3f\n',3000,vp,vs,den);
fclose(fid);
%---------------------------------------------------------------------------------------------


% Generate velocity file for layered model---------------------------------
lmod=dlmread('SMTI-1D-VM.csv',',',1,0); % read in the model parameters
lmod(:,4)=lmod(:,4)*1000; % transfer the unit of density from g/cm^3 to kg/m^3
nl=size(lmod,1); % number of layers
sdep=0; % elevation of the free surface
fid=fopen('layer.dat','wt'); % file name of the velocity model
fprintf(fid,'%9.3f\n',sdep);
for il=1:nl
    % calculate the thickness of each layer
    if il>1
        thickness=lmod(il,1)-lmod(il-1,1);
    else
        thickness=lmod(il,1)-sdep;
    end
    fprintf(fid,'%9.3f  %9.3f  %9.3f %9.3f\n',thickness,lmod(il,2),lmod(il,3),lmod(il,4));
end
fclose(fid);
%--------------------------------------------------------------------------------------------------------------


% Generate velocity model file for 3D modeling-----------------------------
% Note for the original overthrust velocity model, the original storage
% sequence is Z(187)*Y(801)*X(801).
load overthrust_3d_vp.mat; % load 3D data cube of P-wave velocity, m/s

% display the model
n1=1; n2=1; n3=1;
vpdata=permute(data,[2,3,1]); % change the sequence to Y-X-Z
figure; hs=slice(vpdata/1000,n1,n2,n3); hold on;
cp=colorbar('horiz');cp.Label.String='km/s';colormap(jet(512));
axis equal;xlabel('X');ylabel('Y');zlabel('Z');
set(gca,'zdir','reverse'); set(hs,'EdgeColor','none');
set(hs,'DiffuseStrength',0.8); daspect([1,1,1]);
view(-42,16); camproj perspective;

% interpolate the original model according to the model geometry
pmtdata=permute(data,[3 2 1]); % adjust the storage sequence to X-Y-Z
idata=intpmodel3(pmtdata,lx,ly,lz,dx,dy,dz); % note the output data is also in X-Y-Z order

% display the model after inerpolation
n1=1; n2=1; n3=1;
vpdata=permute(idata,[2,1,3]); % change the sequence to Y-X-Z
figure; hs=slice(vpdata/1000,n1,n2,n3); hold on;
cp=colorbar('horiz');cp.Label.String='km/s';colormap(jet(512));
axis equal;xlabel('X');ylabel('Y');zlabel('Z');
set(gca,'zdir','reverse'); set(hs,'EdgeColor','none');
set(hs,'DiffuseStrength',0.8); daspect([1,1,1]);
view(-42,16); camproj perspective;

% generate the S-wave velocity model according to empirical formula
vs=idata/1.7; % let Vp/Vs=1.7, unit: m/s
% generate the density velocity model according to empirical formula
% den=1.741*vp^(0.25)  velocity in km/s, density in g/cm^3 (from Gardner 1974).
den=(1.741*(idata/1000).^0.25)*1000; % now den is in kg/m^3

% output VP binary file for 3D modeling
vp=permute(idata,[3 2 1]); % note the output order should be Z-Y-X according to the 3D modeling
save overthrust_vp.mat vp dx dy dz -v7.3; % save as .mat file
fid1=fopen('VP','w'); % name of the file
fwrite(fid1,vp,'single'); % note the precision, keep consitent with the 3D modeling software
fclose(fid1);

% output VS binary file for 3D modeling
vs=permute(vs,[3 2 1]); % note the output order should be Z-Y-X according to the 3D modeling
save overthrust_vs.mat vs dx dy dz -v7.3; % save as .mat file
fid1=fopen('VS','w'); % name of the file
fwrite(fid1,vs,'single'); % note the precision, keep consitent with the 3D modeling software
fclose(fid1);

% output DEN binary file for 3D modeling
den=permute(den,[3 2 1]); % note the output order should be Z-Y-X according to the 3D modeling
save overthrust_den.mat den dx dy dz -v7.3; % save as .mat file
fid1=fopen('DEN','w'); % name of the file
fwrite(fid1,den,'single'); % note the precision, keep consitent with the 3D modeling software
fclose(fid1);
%--------------------------------------------------------------------------------------------------------------
