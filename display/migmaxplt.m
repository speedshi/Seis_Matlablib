function [tn,xn,yn,zn]=migmaxplt(data,soup,tarxr,taryr,tarzr,taxis)
% This function is used to find the maximum value of the 4D migration data.
% This function also plots the X, Y and Z profiles through the maximum value.
% This function also plots the stacking traces at the located source position.
% This function also plots the X, Y and Z profiles which project maximum
% values along the time dimension.
% The position of the maximum data in the 4D data represents the original
% time and X-Y-Z positions of the source.
% Input:----------------------------------------------------
% data: input 4D migration data. Dimension 1: time; Dimension 2: X;
% Dimension 3: Y; Dimension 4: Z;
% soup: ture source position (X-Y-Z in Km) (dimension: 1*3 or 3*1, a vector);
% tarxr: X range of the target zone (Km) (vector: 1*2 or 2*1);
% taryr: Y range of the target zone (Km) (vector: 1*2 or 2*1);
% tarzr: Z range of the target zone (Km) (vector: 1*2 or 2*1);
% taxis: time axis for showing search origin times (vector: n_times*1);
% Output:--------------------------------------------------
% tn: time index;
% xn: X index;
% yn: Y index;
% zn: Z index.

% set default parameters
if nargin < 6
    taxis = [];
end


% find the position of the maximum value in the 4D data
[~,indx]=max(data(:)); % index of maximum value in the 1D data
[tn,xn,yn,zn]=ind2sub(size(data),indx); % transfer the index to 4D subscribs


% define source location
[~,nxr,nyr,nzr]=size(data);
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
pdis=squeeze(data(tn,xn,:,:))'; % extract the X profile
xx=[taryr(1) taryr(2)];yy=[tarzr(1) tarzr(2)]; % horizontal and vertical range
figure;imagesc(xx,yy,pdis);colormap('jet');colorbar; hold on;
plot(soup(2),soup(3),'kh'); hold on; %caxis([0 1]);
axis equal; axis tight; xlabel('Y (km)'); ylabel('Z (km)');
figure; surf(yyr,zzr,pdis,'LineStyle','none','FaceColor','interp'); hold on;
plot3(soup(2),soup(3),pdis(s1id(3),s1id(2)),'ko','MarkerSize',5,'MarkerFaceColor','k'); hold on;
axis tight; axis off; colormap('jet'); colorbar; title('YZ-plane'); %caxis([0 1]);

% plot the Y profile
pdis=squeeze(data(tn,:,yn,:))'; % extract the Y profile
xx=[tarxr(1) tarxr(2)];yy=[tarzr(1) tarzr(2)]; % horizontal and vertical range
figure;imagesc(xx,yy,pdis);colormap('jet');colorbar; hold on;
plot(soup(1),soup(3),'kh'); hold on; %caxis([0 1]);
axis equal; axis tight; xlabel('X (km)'); ylabel('Z (km)');
figure; surf(xxr,zzr,pdis,'LineStyle','none','FaceColor','interp'); hold on;
plot3(soup(1),soup(3),pdis(s1id(3),s1id(1)),'ko','MarkerSize',5,'MarkerFaceColor','k'); hold on;
axis tight; axis off; colormap('jet'); colorbar; title('XZ-plane'); %caxis([0 1]);

% plot the Z profile
pdis=squeeze(data(tn,:,:,zn));
xx=[taryr(1) taryr(2)];yy=[tarxr(1) tarxr(2)]; % horizontal and vertical range
figure;imagesc(xx,yy,pdis);colormap('jet');colorbar; hold on;
plot(soup(2),soup(1),'kh'); hold on; %caxis([0 1]);
axis equal; axis tight; axis xy; xlabel('Y (km)'); ylabel('X (km)');
figure; surf(yyr,xxr,pdis,'LineStyle','none','FaceColor','interp'); hold on;
plot3(soup(2),soup(1),pdis(s1id(1),s1id(2)),'ko','MarkerSize',5,'MarkerFaceColor','k'); hold on;
axis tight; axis off; colormap('jet'); colorbar; title('XY-plane'); %caxis([0 1]);

% Plot the stacking function at the located source position
figure;
if ~isempty(taxis)
    % has input time axis
    plot(taxis,data(:,xn,yn,zn),'k','linewidth',1.6,'Marker','o','MarkerFaceColor','k','MarkerSize',2);
    xlabel('Time');
else
    % no input time axis
    plot(data(:,xn,yn,zn),'k','linewidth',1.6,'Marker','o','MarkerFaceColor','k','MarkerSize',2);
    xlabel('Time samples');
end
ylabel('Stacked energy'); title('Stacking trace at the located source position');
axis tight;

% Plot the three profiles after projecting maximum values along the time
% dimension.
% projection on Y-Z profile
xx=[taryr(1) taryr(2)];yy=[tarzr(1) tarzr(2)];
imagz=squeeze(max(data,[],1)); % find the maximum value in Time domain
pdis=squeeze(max(imagz,[],1))'; % find the maximum value in X domain
figure;imagesc(xx,yy,pdis);colormap('jet(512)');colorbar; hold on;
plot(soup(2),soup(3),'kh'); hold on;
axis equal; axis tight; xlabel('Y (km)'); ylabel('Z (km)');
figure; surf(yyr,zzr,pdis,'LineStyle','none','FaceColor','interp'); hold on;
plot3(soup(2),soup(3),pdis(s1id(3),s1id(2)),'ko','MarkerSize',5,'MarkerFaceColor','k'); hold on;
axis tight; axis off; colormap('jet'); colorbar; title('YZ-plane'); %caxis([0 1]);

% projection on X-Z profile
xx=[tarxr(1) tarxr(2)];yy=[tarzr(1) tarzr(2)];
imagz=squeeze(max(data,[],1)); % find the maximum value in Time domain
pdis=squeeze(max(imagz,[],2))'; % find the maximum value in Y domain
figure;imagesc(xx,yy,pdis);colormap('jet(512)');colorbar; hold on;
plot(soup(1),soup(3),'kh'); hold on;
axis equal; axis tight; xlabel('X (km)'); ylabel('Z (km)');
figure; surf(xxr,zzr,pdis,'LineStyle','none','FaceColor','interp'); hold on;
plot3(soup(1),soup(3),pdis(s1id(3),s1id(1)),'ko','MarkerSize',5,'MarkerFaceColor','k'); hold on;
axis tight; axis off; colormap('jet'); colorbar; title('XZ-plane'); %caxis([0 1]);

% projection on X-Y profile
xx=[taryr(1) taryr(2)];yy=[tarxr(1) tarxr(2)];
imagz=squeeze(max(data,[],1)); % find the maximum value in Time domain
pdis=squeeze(max(imagz,[],3)); % find the maximum value in Z domain
figure;imagesc(xx,yy,pdis);colormap('jet(512)');colorbar; hold on;
plot(soup(2),soup(1),'kh'); hold on;
axis equal; axis tight; axis xy; xlabel('Y (km)'); ylabel('X (km)');
figure; surf(yyr,xxr,pdis,'LineStyle','none','FaceColor','interp'); hold on;
plot3(soup(2),soup(1),pdis(s1id(1),s1id(2)),'ko','MarkerSize',5,'MarkerFaceColor','k'); hold on;
axis tight; axis off; colormap('jet'); colorbar; title('XY-plane'); %caxis([0 1]);

end