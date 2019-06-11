function mgrsprofdisp(cata,migloc,search,trace,thrsd,bdist)
% This function is used to display the XY, XZ and YZ profiles of the
% location results and also compared with the events in the catalogue.
% Note the unit of the input parameters should be consistent.
% Here X is North direction, Y is East direction, Z is Vertical down.
%
% Unit: meter (m). But display the distance in km!
%
% Input:----------------------------------------------------
% cata: events locations in the catalogue, neca*4, Origin_time-X-Y-Z;
% migloc: our event location results, nevt*5, Origin_time-X-Y-Z-Coherency;
% search: matlab structure, specify the searching (migration) area;
% search.north: vector, 1*2, tarxr, coordinate limit (range) in the X direction;
% search.east: vector, 1*2, taryr, coordinate limit (range) in the Y direction;
% search.depth: vector, 1*2, tarzr, coordinate limit (range) in the Z direction;
% trace: matlab structure, contains station information;
% trace.name: station names;
% trace.north: north coordinates;
% trace.east: east coordinates;
% trace.depth: depth coordinates;
% thrsd: threshold value, only show the points which have migration value
% above the threshold value;
% bdist: boundary value, do not show the points which lie within the limit of boundaries.

% set default values
if nargin==3
    trace=[];
    thrsd=0;
    bdist=0;
elseif nargin==4
    thrsd=0;
    bdist=0;
elseif nargin==5
    bdist=0;
end

if size(migloc,2)<5
    migloc(:,5)=1;
end

% get migration boundary
tarxr=search.north; % X-North range
taryr=search.east; % Y-East range
tarzr=search.depth; % Z-Depth range

% find the located seismic events which fulfill the requirements
migloc=migloc(migloc(:,5)>thrsd & migloc(:,2)>=tarxr(1)+bdist & migloc(:,2)<=tarxr(2)-bdist & migloc(:,3)>=taryr(1)+bdist & migloc(:,3)<=taryr(2)-bdist & migloc(:,4)>=tarzr(1)+bdist & migloc(:,4)<=tarzr(2)-bdist,:);

bsize=16; % the size of the ball in the figure

if ~isempty(trace) && ~isempty(trace.north)
    stname=trace.name; % name of receivers
    recpx=trace.north; % X coordinate of receivers (North)
    recpy=trace.east; % Y coordinate of receivers (East)
    recpz=trace.depth; % Z coordinate of receivers (Down, Depth). Negative value means above the sea leval.
else
    recpx=[];
    recpy=[];
    recpz=[];
end

% X-Y profile
figure;
if ~isempty(recpx)
    plot(recpy/1000,recpx/1000,'kv','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',7); hold on; % plot the stations
    for ir=1:length(stname)
        text(recpy(ir)/1000+3.4,recpx(ir)/1000,stname(ir),'FontWeight','bold','color','k'); % show the names of the stations
    end
end
hold on;
if ~isempty(cata)
    scatter(cata(:,3)/1000,cata(:,2)/1000,bsize,'filled','k'); hold on; % show the catalog events
end
scatter(migloc(:,3)/1000,migloc(:,2)/1000,bsize,migloc(:,1)/3600,'filled'); hold on;
colormap('jet(1024)'); ch=colorbar;
ch.Label.String='Origin time'; ch.Label.FontSize=12;
axis equal tight; box on;
%xlim(taryr/1000); ylim(tarxr/1000);
xlabel('East (km)'); ylabel('North (km)');

% X-Z profile
figure;
if ~isempty(cata)
    scatter(cata(:,2)/1000,cata(:,4)/1000,bsize,'filled','k'); hold on;
end
scatter(migloc(:,2)/1000,migloc(:,4)/1000,bsize,migloc(:,1)/3600,'filled'); hold on;
colormap('jet(1024)'); ch=colorbar;
ch.Label.String='Origin time'; ch.Label.FontSize=12;
set(gca,'ydir','reverse'); set(gca,'color',[0.8 0.8 0.8]);
axis equal; box on;
xlim(tarxr/1000); ylim(tarzr/1000);
xlabel('North (km)'); ylabel('Depth (km)');

% Y-Z profile
figure;
if ~isempty(cata)
    scatter(cata(:,3)/1000,cata(:,4)/1000,bsize,'filled','k'); hold on;
end
scatter(migloc(:,3)/1000,migloc(:,4)/1000,bsize,migloc(:,1)/3600,'filled'); hold on;
colormap('jet(1024)'); ch=colorbar;
ch.Label.String='Origin time'; ch.Label.FontSize=12;
set(gca,'ydir','reverse'); set(gca,'color',[0.8 0.8 0.8]);
axis equal; box on;
xlim(taryr/1000); ylim(tarzr/1000);
xlabel('East (km)'); ylabel('Depth (km)');

end