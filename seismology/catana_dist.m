function catana_dist(file,earthquake)
% This function is used to analyse the catalog.
% Show the distances of events in the catalog with respect to the input earthquake.
%
% INPUT--------------------------------------------------------------------
% file: string, file name of the catalog;
% earthquake: structure, information of the input earthquake;
% earthquake.latitude: latitude in degree of the input earthquake;
% earthquake.longitude: longitude in degree of the input earthquake;
% earthquake.elevation: elevation in meter of the input earthquake;
%
% OUTPUT-------------------------------------------------------------------
% figure showing the distance.


% obtain the coordinates of the input earthquake
[east,north,depth]=geod2cart(earthquake.latitude,earthquake.longitude,earthquake.elevation);

% read in the catalog data
catalog=read_catalog(file);


% calculate the distances between the events and the input earthquake
dist=sqrt((catalog.east-east).^2+(catalog.north-north).^2)/1000; % note the unit is now in km

sz=dnormlz(catalog.magnitude,4,12); % magnitude of events control the relative sizes of the ploted dots

% plot the distances
figure;
scatter(catalog.time,dist,sz,catalog.depth/1000,'filled');
cb=colorbar; cb.Label.String='Event depth (km)';
axis tight;
xlabel('Time');
ylabel('Distance from epicenter (km)');
    

end