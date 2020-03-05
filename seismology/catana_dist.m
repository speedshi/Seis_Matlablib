function catana_dist(catalog,earthquake,para)
% This function is used to analyse the catalog.
% Show the distances of events in the catalog with respect to the input earthquake.
%
% INPUT--------------------------------------------------------------------
% catalog: structure, the catalog of seismic events, unit: meter and degree;
% earthquake: structure, information of the input earthquake;
% earthquake.latitude: latitude in degree of the input earthquake;
% earthquake.longitude: longitude in degree of the input earthquake;
% earthquake.elevation: elevation in meter of the input earthquake;
% para: structure, contain parameters to adjust the plot;
% para.depth: if use color to respresent event depth;
% para.newfig: if plot on a new figure;
%
% OUTPUT-------------------------------------------------------------------
% figure showing the distance.

% set default parameters
if nargin<3
    para.depth=true;
    para.newfig=true;
end

if ~isfield(para,'depth')
    para.depth=true;
end

if ~isfield(para,'newfig')
    para.newfig=true;
end


% obtain the coordinates of the input earthquake
[east,north,depth]=geod2cart(earthquake.latitude,earthquake.longitude,earthquake.elevation);

% calculate the distances between the events and the input earthquake
dist=sqrt((catalog.east-east).^2+(catalog.north-north).^2)/1000; % note the unit is now in km

sz=dnormlz(catalog.magnitude,5,20); % magnitude of events control the relative sizes of the ploted dots

% plot the distances
if para.newfig
    figure;
end

hold on;

if para.depth
    scatter(catalog.time,dist,sz,catalog.depth/1000,'filled','MarkerFaceAlpha',0.4,'MarkerEdgeAlpha',0.4);
    cb=colorbar; cb.Label.String='Event depth (km)';
else
    scatter(catalog.time,dist,sz,'k','filled','MarkerFaceAlpha',0.4,'MarkerEdgeAlpha',0.4);
end
axis tight;
ax=gca;
ax.FontSize=14;
xlabel('Time','fontsize',16);
xtickformat('yyyy-MM');
ax.XAxis.MinorTick = 'on';
set(ax,'TickLength',[0.012 0.02]);
ax.LineWidth=1.3;
ylabel('Distance from epicenter (km)','fontsize',16);


end