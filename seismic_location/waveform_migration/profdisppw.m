function profdisppw(cata,migloc,tarxr,taryr,tarzr,eid)
% This function is used to display the XY, XZ and YZ profiles of the
% location results and also compared with the events in the catalogue.
% The events are displayed and compared pair-wisely. So this required the
% 'cata' and 'migloc' have the same number of seismic events that can be
% matched accordingly.
% Note the unit of the input parameters should be consistent.
% Here X is North direction, Y is East direction, Z is Vertical down.
% Input:----------------------------------------------------
% cata: events locations in the catalogue, neca*4, Origin_time-X-Y-Z;
% migloc: our event location results, neca*5, Origin_time-X-Y-Z-Coherency;
% tarxr: coordinate limit (range) in the X direction;
% taryr: coordinate limit (range) in the Y direction;
% tarzr: coordinate limit (range) in the Z direction;
% eid: event index which specifies which events need to be displayed,
% a vector contains the index of events, index must smaller than 'neca'.

neca=size(cata,1);

% set default value
if nargin<=5
    eid=1:1:neca;
end

ne=length(eid); % number of seismic events

% X-Y profile
figure;
for ii=1:ne
    scatter(cata(eid(ii),3),cata(eid(ii),2),'filled','k'); hold on;
    text(cata(eid(ii),3),cata(eid(ii),2),sprintf('%d',eid(ii))); hold on;
    scatter(migloc(eid(ii),3),migloc(eid(ii),2),36,migloc(eid(ii),5),'filled'); hold on; 
    colormap(jet);colorbar;
    text(migloc(eid(ii),3),migloc(eid(ii),2),sprintf('%d',eid(ii))); hold on;
end
axis equal; xlim(taryr); ylim(tarxr); box on;
xlabel('East (km)'); ylabel('North (km)');

% X-Z profile
figure;
for ii=1:ne
    scatter(cata(eid(ii),2),cata(eid(ii),4),'filled','k'); hold on;
    text(cata(eid(ii),2),cata(eid(ii),4),sprintf('%d',eid(ii))); hold on;
    scatter(migloc(eid(ii),2),migloc(eid(ii),4),'filled','b'); hold on;
    text(migloc(eid(ii),2),migloc(eid(ii),4),sprintf('%d',eid(ii))); hold on;
end
set(gca,'Ydir','reverse');
axis equal; xlim(tarxr); ylim(tarzr); box on;
xlabel('North (km)'); ylabel('Depth (km)');

% Y-Z profile
figure;
for ii=1:ne
    scatter(cata(eid(ii),3),cata(eid(ii),4),'filled','k'); hold on;
    text(cata(eid(ii),3),cata(eid(ii),4),sprintf('%d',eid(ii))); hold on;
    scatter(migloc(eid(ii),3),migloc(eid(ii),4),'filled','b'); hold on;
    text(migloc(eid(ii),3),migloc(eid(ii),4),sprintf('%d',eid(ii))); hold on;
end
set(gca,'Ydir','reverse');
axis equal; xlim(taryr); ylim(tarzr); box on;
xlabel('East (km)'); ylabel('Depth (km)');

end