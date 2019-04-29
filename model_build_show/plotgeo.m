function plotgeo(xc,yc,zc,lydp,recx,recy,recz,rwlxy,rwlz,srar,soup,rmksize,smksize,wmksize)
% This function is used to plot the geometry of the model and recording arrays.
% Right-hand coordinate systerm in Cartesian coordinate, with X-North,
% Y-East, Z-Vertical down. Unit: meter.
% INPUT--------------------------------------------------------
% xc: boundry of X axis of the model (m), vector (2*1), can be void vector;
% yc: boundry of Y axis of the model (m), vector (2*1), can be void vector;
% zc: boundry of Z axis of the model (m), vector (2*1), can be void vector;
% lydp: depth of each layer interface, should include the inerface of the 
% free surface (m). For free surface: 0-->sea-level, (+)positive-->below
% sea-level, (-)negative-->above sea-level;
% recx: X coordinates of surface arrays (m), a vector, nr*1;
% recy: Y coordinates of surface arrays (m), a vector, nr*1;
% recz: Z coordinates of surface arrays (m), a vector, nr*1;
% rwlxy: hrizontal position of downhole arrays(X-Y) (m), nwell*2, could be
% more than one downhole array, then each row is an array position; note
% the downhole array may not exist, thus 'rwlxy' can be a null matrix;
% rwlz: Z positon range of downhole arrays (m), must has the same rows as 'rwlxy';
% sara: shows the migration searching area (m), is 3*2 matrix:
% row 1: X range; row 2: Y range; row 3: Z range;
% soup: position of source points(X-Y-Z) (m), ns*3, could be more than
% one source, each row is a source;
% rmksize: markersize of the surface receivers;
% smksize: markersize of the sources;
% wmksize: markersize of the downhole receivers.
% OUTPUT---------------------------------------------------------
% Three figures show the model geometry and surface projection
% and 2 vertical projection of sources and downhole arrays.

if nargin<8
    rwlxy=[];
    rwlz=[];
    srar=[];
    soup=[];
    rmksize=6;
    smksize=6;
    wmksize=4;
elseif nargin<10
    srar=[];
    soup=[];
    rmksize=6;
    smksize=6;
    wmksize=4;
elseif nargin<11
    soup=[];
    rmksize=6;
    smksize=6;
    wmksize=4;
elseif nargin<12
    rmksize=6;
    smksize=6;
    wmksize=4;
end

if isempty(rwlxy)
    wxmin=[];
    wxmax=[];
    wymin=[];
    wymax=[];
else
    wxmin=min(rwlxy(:,1));
    wxmax=max(rwlxy(:,1));
    wymin=min(rwlxy(:,2));
    wymax=max(rwlxy(:,2));
end

if isempty(srar)
    axmin=[];
    axmax=[];
    aymin=[];
    aymax=[];
    azmin=[];
    azmax=[];
else
    axmin=min(srar(1,:));
    axmax=max(srar(1,:));
    aymin=min(srar(2,:));
    aymax=max(srar(2,:));
    azmin=min(srar(3,:));
    azmax=max(srar(3,:));
end

if isempty(soup)
    sxmin=[];
    sxmax=[];
    symin=[];
    symax=[];
    szmin=[];
    szmax=[];
else
    sxmin=min(soup(:,1));
    sxmax=max(soup(:,1));
    symin=min(soup(:,2));
    symax=max(soup(:,2));
    szmin=min(soup(:,3));
    szmax=max(soup(:,3));
end

if isempty(xc)
    min_temp=min([recx(:); wxmin; sxmin; axmin]);
    max_temp=max([recx(:); wxmax; sxmax; axmax]);
    xc(1)= min_temp-0.1*abs(max_temp-min_temp);    
    xc(2)= max_temp+0.1*abs(max_temp-min_temp);
end

if isempty(yc)
    min_temp=min([recy(:); wymin; symin; aymin]);
    max_temp=max([recy(:); wymax; symax; aymax]);
    yc(1)= min_temp-0.1*abs(max_temp-min_temp);    
    yc(2)= max_temp+0.1*abs(max_temp-min_temp);
end

if isempty(zc)
    min_temp=min([recz(:); min(rwlz(:)); szmin; azmin; lydp(1)]);
    max_temp=max([recz(:); max(rwlz(:)); szmax; azmax; lydp(end)]);
    zc(1)= min_temp;    
    zc(2)= max_temp+0.1*abs(max_temp-min_temp);
end

% transfer meters to kilometers
xc=xc/1000;
yc=yc/1000;
zc=zc/1000;
lydp=lydp/1000;
recx=recx/1000;
recy=recy/1000;
recz=recz/1000;
rwlxy=rwlxy/1000;
rwlz=rwlz/1000;
srar=srar/1000;
soup=soup/1000;

if ~isempty(rwlxy)
    rwlpx=rwlxy(:,1); % X positon of downhole arrays
    rwlpy=rwlxy(:,2); % Y positon of downhole arrays
end
ns=size(soup,1); % number of sources
ndh=size(rwlxy,1); % number of downhole arrays

xg=[xc(1) xc(1) xc(2) xc(2) xc(1)];
yg=[yc(1) yc(2) yc(2) yc(1) yc(1)];
zg1=zc(1)*ones(size(xg));
zg2=zc(2)*ones(size(xg));

lyxl=[xc(1) xc(1) xc(2)];
lyyl=[yc(1) yc(2) yc(2)];
lyxd=[xc(2) xc(2) xc(1)];
lyyd=[yc(2) yc(1) yc(1)];
nly=max(size(lydp));

figure;
plot3(xg,yg,zg1,'k',xg,yg,zg2,'k');hold on; % plot the model
axis tight;
if ~isempty(soup)
    for is=1:ns
        plot3(soup(is,1),soup(is,2),soup(is,3),'rp','markersize',smksize);hold on; % plot the source point
    end
end
plot3(recx,recy,recz,'b.','markersize',rmksize);hold on; % plot the surface arrays
if ~isempty(rwlxy)
    for iw=1:ndh
        rwlx=rwlpx(iw)*ones(size(rwlz(iw,:)));
        rwly=rwlpy(iw)*ones(size(rwlz(iw,:)));
        plot3(rwlx,rwly,rwlz(iw,:),'bv','markersize',wmksize);hold on; % plot the down hole arrays
    end
end
for i=1:nly
    % plot the layers
    lyz=lydp(i)*ones(size(lyxl));
    plot3(lyxl,lyyl,lyz,'k');hold on;
    plot3(lyxd,lyyd,lyz,':k');hold on;
end
xlabel('X (km)');ylabel('Y (km)');zlabel('Z (km)');
box on;set(gca,'BoxStyle','full');
set(gca,'Xdir','reverse');set(gca,'Zdir','reverse');
set(gca,'XMinorGrid','on');set(gca,'YMinorGrid','on');set(gca,'ZMinorGrid','on');
set(gca,'MinorGridLineStyle',':');
if ~isempty(srar)
    cxx=[srar(1,1) srar(1,1) srar(1,1) srar(1,2) srar(1,2) srar(1,1);srar(1,1) srar(1,2) srar(1,1) srar(1,1) srar(1,2) srar(1,1);srar(1,1) srar(1,2) srar(1,2) srar(1,1) srar(1,2) srar(1,2);srar(1,1) srar(1,1) srar(1,2) srar(1,2) srar(1,2) srar(1,2)];
    cyy=[srar(2,1) srar(2,2) srar(2,1) srar(2,1) srar(2,1) srar(2,1);srar(2,2) srar(2,2) srar(2,2) srar(2,1) srar(2,2) srar(2,2);srar(2,2) srar(2,2) srar(2,2) srar(2,1) srar(2,2) srar(2,2);srar(2,1) srar(2,2) srar(2,1) srar(2,1) srar(2,1) srar(2,1)];
    czz=[srar(3,1) srar(3,1) srar(3,1) srar(3,1) srar(3,1) srar(3,2);srar(3,1) srar(3,1) srar(3,1) srar(3,1) srar(3,1) srar(3,2);srar(3,2) srar(3,2) srar(3,1) srar(3,2) srar(3,2) srar(3,2);srar(3,2) srar(3,2) srar(3,1) srar(3,2) srar(3,2) srar(3,2)];
    fill3(cxx,cyy,czz,'y','linestyle','none','FaceAlpha',0.35); hold on;
end
view([116,22]);

% plot the surface projection of sources, surface and downhole arrays
figure;
plot(recx,recy,'b.','markersize',rmksize); hold on; % plot the projection of surface arrays
if ~isempty(soup)
    for is=1:ns
        plot(soup(is,1),soup(is,2),'rp','markersize',smksize);hold on; % sources
    end
end
if ~isempty(rwlxy)
    for iw=1:ndh
        plot(rwlpx(iw),rwlpy(iw),'bv','markersize',wmksize);hold on; % downhole arrays
    end
end
axis equal;
xlim([xc(1) xc(2)]);
ylim([yc(1) yc(2)]);
if ~isempty(srar)
    fill([srar(1,1) srar(1,2) srar(1,2) srar(1,1)],[srar(2,1) srar(2,1) srar(2,2) srar(2,2)],'y','linestyle','none','FaceAlpha',0.35);
end
axis ij;xlabel('X (km)');ylabel('Y (km)');
set(gca,'XMinorGrid','on');set(gca,'YMinorGrid','on');
view(-90,90); % vertical direction is X axis, horizontal direction is Y axis.

% plot the vertical projection of sources, surface and downhole arrays
% project on the middle of the model and perpendicular to X axis.
figure;
plot(recy,recz,'b.','markersize',rmksize);hold on; % plot the projection of surface arrays
if ~isempty(soup)
    for is=1:ns
        plot(soup(is,2),soup(is,3),'rp','markersize',smksize);hold on; % sources
    end
end
if ~isempty(rwlxy)
    for iw=1:ndh
        rwly=rwlpy(iw)*ones(size(rwlz(iw,:)));
        plot(rwly,rwlz(iw,:),'bv','markersize',wmksize);hold on; % downhole arrays
    end
end
for i=1:nly
    plot([yc(1) yc(2)],[lydp(i) lydp(i)],'k','linewidth',1.6);hold on;% plot the layer boundarys
end
axis equal;
xlim([yc(1) yc(2)]); ylim([zc(1) zc(2)]);
if ~isempty(srar)
    fill([srar(2,1) srar(2,2) srar(2,2) srar(2,1)],[srar(3,1) srar(3,1) srar(3,2) srar(3,2)],'y','linestyle','none','FaceAlpha',0.35);
end
axis ij;xlabel('Y (km)');ylabel('Z (km)');
set(gca,'XMinorGrid','on');set(gca,'YMinorGrid','on');

% plot the vertical projection of sources, surface and downhole arrays
% project on the middle of the model and perpendicular to Y axis.
figure;
plot(recx,recz,'b.','markersize',rmksize);hold on; % plot the projection of surface arrays
if ~isempty(soup)
    for is=1:ns
        plot(soup(is,1),soup(is,3),'rp','markersize',smksize);hold on; % sources
    end
end
if ~isempty(rwlxy)
    for iw=1:ndh
        rwlx=rwlpx(iw)*ones(size(rwlz(iw,:)));
        plot(rwlx,rwlz(iw,:),'bv','markersize',wmksize);hold on; % downhole arrays
    end
end
for i=1:nly
    plot([xc(1) xc(2)],[lydp(i) lydp(i)],'k','linewidth',1.6);hold on;% plot the layer boundarys
end
axis equal;
xlim([xc(1) xc(2)]); ylim([zc(1) zc(2)]);
if ~isempty(srar)
    fill([srar(1,1) srar(1,2) srar(1,2) srar(1,1)],[srar(3,1) srar(3,1) srar(3,2) srar(3,2)],'y','linestyle','none','FaceAlpha',0.35);
end
axis ij;xlabel('X (km)');ylabel('Z (km)');
set(gca,'XMinorGrid','on');set(gca,'YMinorGrid','on');

end

