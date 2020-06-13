function [tn,xn,yn,zn]=migmaxplt(data,soup,tarxr,taryr,tarzr,para)
% This function is used to find the maximum value of the 4D migration data.
% This function also plots the X, Y and Z profiles through the maximum value.
% This function also plots the stacking traces at the located source position.
% This function also plots the X, Y and Z profiles which project maximum
% values along the time dimension.
% The position of the maximum data in the 4D data represents the original
% time and X-Y-Z positions of the source.
% X-North; Y-East; Z-Depth;
% Input:----------------------------------------------------
% data: input 4D migration data. Dimension 1: time; Dimension 2: X;
% Dimension 3: Y; Dimension 4: Z;
% soup: ture source position (X-Y-Z in Km) (dimension: 1*3 or 3*1, a vector);
% tarxr: X range of the target zone (Km) (vector: 1*2 or 2*1);
% taryr: Y range of the target zone (Km) (vector: 1*2 or 2*1);
% tarzr: Z range of the target zone (Km) (vector: 1*2 or 2*1);
% para: structure, controlling plotting parameters;
% para.taxis: vector: n_times*1, time axis for showing search origin times;
% para.ctlpct: scalar, plot a contour-line on the migration profile;
% Output:------------------------------------------------------------------
% tn: time index;
% xn: X index;
% yn: Y index;
% zn: Z index.


% set default parameters---------------------------------------------------
if nargin < 6
    para.taxis = [];
    para.ctlpct = [];
end

if ~isfield(para,'taxis')
    para.taxis = [];
end

if ~isfield(para,'ctlpct')
    para.ctlpct = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% find the position of the maximum value in the 4D data
[mgvmax,indx]=max(data(:)); % index of maximum value in the 1D data
[tn,xn,yn,zn]=ind2sub(size(data),indx); % transfer the index to 4D subscribs

if ~isempty(para.ctlpct)
    ctlv = para.ctlpct*mgvmax; % determine the contour-line level
end

% define source location
[nt0,nxr,nyr,nzr]=size(data);
xxr=linspace(tarxr(1),tarxr(2),nxr);
yyr=linspace(taryr(1),taryr(2),nyr);
zzr=linspace(tarzr(1),tarzr(2),nzr);
if ~isempty(soup)
    % have input for source location
    [~,s1id(1)]=min(abs(xxr-soup(1)));
    [~,s1id(2)]=min(abs(yyr-soup(2)));
    [~,s1id(3)]=min(abs(zzr-soup(3)));
else
    % no input for source location, define the source using the maximum
    % migration value in the data
    s1id=[xn,yn,zn]; % the index of the source, i.e. the point having the maximum migration value
    soup=[xxr(xn),yyr(yn),zzr(zn)]; % the location of the source
end


% load a specified colormap
load mycolor1.mat;

% Plot the three profiles along the maximum migration value in the volume
% plot the X profile
pdis=reshape(data(tn,xn,:,:),nyr,nzr)'; % extract the X profile
if ~isvector(pdis)
    xx=[taryr(1) taryr(2)];yy=[tarzr(1) tarzr(2)]; % horizontal and vertical range
    figure;imagesc(xx,yy,pdis);colormap('jet');colorbar; hold on;
    plot(soup(2),soup(3),'kp','MarkerSize',9); hold on; %caxis([0 1]);
    if ~isempty(para.ctlpct)
        cxx = linspace(xx(1),xx(2),size(pdis,1));
        cyy = linspace(yy(1),yy(2),size(pdis,2));
        contour(gca,cxx,cyy,pdis,[ctlv ctlv],'k','linewidth',1.2); % plot the contour-line
    end
    axis equal; axis tight; xlabel('East (km)'); ylabel('Depth (km)'); title('Sliced profile');
    
    figure; surf(yyr,zzr,pdis,'LineStyle','none','FaceColor','interp'); hold on;
    plot3(soup(2),soup(3),pdis(s1id(3),s1id(2)),'ko','MarkerSize',6,'MarkerFaceColor','k'); hold on;
    axis tight; axis off; colormap('jet'); colorbar; title('East-Depth plane (sliced)'); %caxis([0 1]);
end

% plot the Y profile
pdis=reshape(data(tn,:,yn,:),nxr,nzr)'; % extract the Y profile
if ~isvector(pdis)
    xx=[tarxr(1) tarxr(2)];yy=[tarzr(1) tarzr(2)]; % horizontal and vertical range
    figure;imagesc(xx,yy,pdis);colormap('jet');colorbar; hold on;
    plot(soup(1),soup(3),'kp','MarkerSize',9); hold on; %caxis([0 1]);
    axis equal; axis tight; xlabel('North (km)'); ylabel('Depth (km)'); title('Sliced profile');
    if ~isempty(para.ctlpct)
        cxx = linspace(xx(1),xx(2),size(pdis,1));
        cyy = linspace(yy(1),yy(2),size(pdis,2));
        contour(gca,cxx,cyy,pdis,[ctlv ctlv],'k','linewidth',1.2); % plot the contour-line
    end
    
    figure; surf(xxr,zzr,pdis,'LineStyle','none','FaceColor','interp'); hold on;
    plot3(soup(1),soup(3),pdis(s1id(3),s1id(1)),'ko','MarkerSize',6,'MarkerFaceColor','k'); hold on;
    axis tight; axis off; colormap('jet'); colorbar; title('North-Depth plane (sliced)'); %caxis([0 1]);
end

% plot the Z profile
pdis=reshape(data(tn,:,:,zn),nxr,nyr);
if ~isvector(pdis)
    xx=[taryr(1) taryr(2)];yy=[tarxr(1) tarxr(2)]; % horizontal and vertical range
    figure;imagesc(xx,yy,pdis);colormap('jet');colorbar; hold on;
    plot(soup(2),soup(1),'kp','MarkerSize',9); hold on; %caxis([0 1]);
    axis equal; axis tight; axis xy; xlabel('East (km)'); ylabel('North (km)'); title('Sliced profile');
    if ~isempty(para.ctlpct)
        cxx = linspace(xx(1),xx(2),size(pdis,1));
        cyy = linspace(yy(1),yy(2),size(pdis,2));
        contour(gca,cxx,cyy,pdis,[ctlv ctlv],'k','linewidth',1.2); % plot the contour-line
    end
    
    figure; surf(yyr,xxr,pdis,'LineStyle','none','FaceColor','interp'); hold on;
    plot3(soup(2),soup(1),pdis(s1id(1),s1id(2)),'ko','MarkerSize',6,'MarkerFaceColor','k'); hold on;
    axis tight; axis off; colormap('jet'); colorbar; title('East-North plane (sliced)'); %caxis([0 1]);
end

% Plot the stacking function at the located source position
figure;
if ~isempty(para.taxis)
    % has input time axis
    plot(para.taxis,data(:,xn,yn,zn),'k','linewidth',1.6,'Marker','o','MarkerFaceColor','k','MarkerSize',2);
    hold on;
    plot(para.taxis(tn),data(tn,xn,yn,zn),'r','Marker','o','MarkerFaceColor','r','MarkerSize',6);
    xlabel('Time');
else
    % no input time axis
    plot(data(:,xn,yn,zn),'k','linewidth',1.6,'Marker','o','MarkerFaceColor','k','MarkerSize',2); 
    hold on;
    plot(tn,data(tn,xn,yn,zn),'r','Marker','o','MarkerFaceColor','r','MarkerSize',6);
    xlabel('Time samples');
end
ylabel('Stacked energy'); title('Stacking trace at the located source position');
axis tight;

% Plot the three profiles after projecting maximum values along the time
% dimension.
% projection on Y-Z profile
xx=[taryr(1) taryr(2)];yy=[tarzr(1) tarzr(2)];
imagz=reshape(max(data,[],1),nxr,nyr,nzr); % find the maximum value in Time domain
pdis=reshape(max(imagz,[],1),nyr,nzr)'; % find the maximum value in X domain
if ~isvector(pdis)
    figure;imagesc(xx,yy,pdis);colormap('jet(512)');colorbar; hold on;
    plot(soup(2),soup(3),'kp','MarkerSize',9); hold on;
    axis equal; axis tight; xlabel('East (km)'); ylabel('Depth (km)'); title('Projected profile');
    if ~isempty(para.ctlpct)
        cxx = linspace(xx(1),xx(2),size(pdis,1));
        cyy = linspace(yy(1),yy(2),size(pdis,2));
        contour(gca,cxx,cyy,pdis,[ctlv ctlv],'k','linewidth',1.2); % plot the contour-line
    end
    
    figure; surf(yyr,zzr,pdis,'LineStyle','none','FaceColor','interp'); hold on;
    plot3(soup(2),soup(3),pdis(s1id(3),s1id(2)),'ko','MarkerSize',6,'MarkerFaceColor','k'); hold on;
    axis tight; axis off; colormap('jet'); colorbar; title('East-Depth plane (projected)'); %caxis([0 1]);
end

% projection on X-Z profile
xx=[tarxr(1) tarxr(2)];yy=[tarzr(1) tarzr(2)];
imagz=reshape(max(data,[],1),nxr,nyr,nzr); % find the maximum value in Time domain
pdis=reshape(max(imagz,[],2),nxr,nzr)'; % find the maximum value in Y domain
if ~isvector(pdis)
    figure;imagesc(xx,yy,pdis);colormap('jet(512)');colorbar; hold on;
    plot(soup(1),soup(3),'kp','MarkerSize',9); hold on;
    axis equal; axis tight; xlabel('North (km)'); ylabel('Depth (km)'); title('Projected profile');
    if ~isempty(para.ctlpct)
        cxx = linspace(xx(1),xx(2),size(pdis,1));
        cyy = linspace(yy(1),yy(2),size(pdis,2));
        contour(gca,cxx,cyy,pdis,[ctlv ctlv],'k','linewidth',1.2); % plot the contour-line
    end
    
    figure; surf(xxr,zzr,pdis,'LineStyle','none','FaceColor','interp'); hold on;
    plot3(soup(1),soup(3),pdis(s1id(3),s1id(1)),'ko','MarkerSize',6,'MarkerFaceColor','k'); hold on;
    axis tight; axis off; colormap('jet'); colorbar; title('North-Depth plane (projected)'); %caxis([0 1]);
end

% projection on X-Y profile
xx=[taryr(1) taryr(2)];yy=[tarxr(1) tarxr(2)];
imagz=reshape(max(data,[],1),nxr,nyr,nzr); % find the maximum value in Time domain
pdis=reshape(max(imagz,[],3),nxr,nyr); % find the maximum value in Z domain
if ~isvector(pdis)
    figure;imagesc(xx,yy,pdis);colormap('jet(512)');colorbar; hold on;
    plot(soup(2),soup(1),'kp','MarkerSize',9); hold on;
    axis equal; axis tight; axis xy; xlabel('East (km)'); ylabel('North (km)'); title('Projected profile');
    if ~isempty(para.ctlpct)
        cxx = linspace(xx(1),xx(2),size(pdis,1));
        cyy = linspace(yy(1),yy(2),size(pdis,2));
        contour(gca,cxx,cyy,pdis,[ctlv ctlv],'k','linewidth',1.2); % plot the contour-line
    end
    
    figure; surf(yyr,xxr,pdis,'LineStyle','none','FaceColor','interp'); hold on;
    plot3(soup(2),soup(1),pdis(s1id(1),s1id(2)),'ko','MarkerSize',6,'MarkerFaceColor','k'); hold on;
    axis tight; axis off; colormap('jet'); colorbar; title('East-North plane (projected)'); %caxis([0 1]);
end

end