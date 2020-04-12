function plotmodel(sname,rname,mname)
% This function shows the model.
% INPUT:-------------------------------------------------
% sname: the file name of source file;
% rname: the file name of receiver file;
% mname: the file name of velocity model file.

if nargin < 1
    sname="source.dat";
    rname="receiver.dat";
    mname="model.dat";
end

% read the source file
[~, ~, ns, source]=rdsourcef(sname);
soup=zeros(ns,3); % obtain source locations
for ids=1:ns
    soup(ids,:)=source(ids).pos;
end

% read the receiver file
stations=rdreceiverf(rname);
recp=stations.recp; % positions of stations, X-Y-Z

% read the velocity model file
[nl,rsd0,model]=rdmodelf(mname);
% calculate the depth of each internal layer interface, exclude the top layer (free surface) and the bottom layer.
if nl==1
    % only one layer (homogeneous model), no interval layer interface
    lydp=[];
elseif nl==2
    lydp(1)=rsd0+model.thickness(1);
else
    lydp=zeros(nl-1,1);
    lydp(1)=rsd0+model.thickness(1);
    for id=2:nl-1
        lydp(id)=lydp(id-1)+model.thickness(id);
    end
end

% determine the boundry of the model
xmin=min([min(recp(:,1)) min(soup(:,1))]);
xmax=max([max(recp(:,1)) max(soup(:,1))]);
ymin=min([min(recp(:,2)) min(soup(:,2))]);
ymax=max([max(recp(:,2)) max(soup(:,2))]);
zmin=min([min(recp(:,3)) min(soup(:,3)) rsd0]);
zmax=max([max(recp(:,3)) max(soup(:,3)) rsd0+sum(model.thickness(:))]);

xxc=0.2*(xmax-xmin);
xc=[xmin-xxc xmax+xxc];
yyc=0.2*(ymax-ymin);
yc=[ymin-yyc ymax+yyc];
zzc=0.05*(zmax-zmin);
zc=[zmin zmax+zzc];

para.sname = stations.name;
plotgeo(xc,yc,zc,lydp,recp(:,1),recp(:,2),recp(:,3),[],[],[],soup,16,10,8,para);

end