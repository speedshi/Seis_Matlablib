function plot_evesta(stations,events)
% This function is used to display the stations and seismic events.
% Note if you input a matlab structure, it is better not to modifile it in
% the function. Because that will modified the structure event outside the
% function.
%
% INPUT--------------------------------------------------------------------
% stations: structure, contains the station information;
% stations.name: cell array, the name of each station;
% stations.latitude: vector, latitude in degree of each station;
% stations.longitude: vector, longitude in degree of each station;
% events: structure, contains the seismic event information;
% events.latitude: vector, latitude in degree of each event;
% events.longitude: vector, longitude in degree of each event;


% obtain Cartisian coordinates of stations
[sta.east,sta.north,~]=geod2cart(stations.latitude,stations.longitude,0);

% obtain Cartisian coordinates of seismic events
[eve.east,eve.north,~]=geod2cart(events.latitude,events.longitude,0);

east_max=max([sta.east(:)/1000; eve.east(:)/1000]);
east_min=min([sta.east(:)/1000; eve.east(:)/1000]);
east_d=0.16*(east_max-east_min);

north_max=max([sta.north(:)/1000; eve.north(:)/1000]);
north_min=min([sta.north(:)/1000; eve.north(:)/1000]);
north_d=0.16*(north_max-north_min);

size_sta=28; % marker size of the stations
size_eve=28; % marker size of the events

figure;
scatter(sta.east/1000,sta.north/1000,size_sta,'b^','filled'); hold on; % plot stations
scatter(eve.east/1000,eve.north/1000,size_eve,'rp','filled'); hold on; % plot events
text(sta.east/1000+2,sta.north/1000,stations.name,'FontSize',10); % show station name
axis equal; box on;
xlim([east_min-east_d east_max+east_d]);
ylim([north_min-north_d north_max+north_d]);
xlabel('East (km)'); ylabel('North (km)');


end