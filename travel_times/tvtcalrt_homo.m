function [tvtp,tvts,tkfa,rayl]=tvtcalrt_homo(vp,vs,soup,recp)
% This function is used to calculate the travel-times, take-off angles and
% ray-length in homogeneous isotropic media. The source is placed at target
% layer, receiver is placed at the free surface or near free surface. The
% receiver position locate above the source position.
% Usage:
% input:----------------------------------------------------
% 1 vp: velocity of P-wave, unit: km/s;
% 2 vs: velocity of S-wave, unit: km/s;
% 3 soup: source position, dimension: nsr*3, each row is a source position,
% column 1: X positon; column 2: Y position; column 3: Z position; unit: km;
% 4 recp: receiver position, dimension: nre*3,each row is a receiver position,
% column 1: X positon; column 2: Y position; column 3: Z position, unit: km.
% output:--------------------------------------------------
% 1 tvtp: travel-time of direct P-wave, dimension: nsr*nre, row corresponds to source index,
% column corresponds to receiver index; unit: second.
% 2 tvts: travel-time of direct S-wave, dimension: nsr*nre, row corresponds to source index,
% column corresponds to receiver index; unit: second.
% 3 tkfa: take-off angle, dimension: nsr*nre, row corresponds to source index,
% column corresponds to receiver index; unit: degree; Take-off angle: the
% angle of rays with respect to vertical direction, range: [0 90]; equal to
% incident angle of the rays here;
% 4 rayl: length of the ray, dimension: nsr*nre, row corresponds to source index,
% column corresponds to receiver index, unit: km.

nsr=size(soup,1); % number of sources
nre=size(recp,1); % number of surface geophones

% initialization
tvtp=zeros(nsr,nre);
tvts=zeros(nsr,nre);
tkfa=zeros(nsr,nre);
rayl=zeros(nsr,nre);

recp_x = recp(:,1);  % X-coordinates of stations
recp_y = recp(:,2);  % Y-coordinates of stations
recp_z = recp(:,3);  % Z-coordinates of stations

soup_x = soup(:,1);  % X-coordinates of sources
soup_y = soup(:,2);  % Y-coordinates of sources
soup_z = soup(:,3);  % Z-coordinates of sources

% calculate the ray-length of the source-receiver pairs and travel-times
parfor ir=1:nre
    for is=1:nsr
        % calculate the source-receiver offset
        srx=recp_x(ir)-soup_x(is); % source-receiver distance in X direction (km)
        sry=recp_y(ir)-soup_y(is); % source-receiver distance in Y direction (km)
        srz=recp_z(ir)-soup_z(is); % source-receiver distance in Z direction (km)
        sroff=sqrt(srx*srx+sry*sry); % horizontal distance of source-receiver (km)
        srdist=sqrt(srx*srx+sry*sry+srz*srz); % source-receiver distance (km)        
        tvtp(is,ir)=srdist/vp; % travel-time of direct P-wave
        tvts(is,ir)=srdist/vs; % travel-time of direct S-wave
        rayl(is,ir)=srdist; % ray-length
        tkfa(is,ir)=asind(sroff/srdist); % take-off angles       
    end
end

end