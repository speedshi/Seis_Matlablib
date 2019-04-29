function [tvt,tkfa,rayl]=tvtcalrt_ly(vel,thk,soup,recp)
% This function is used to calculate the travel-times, take-off angles and 
% ray-length in layered isotropic media using ray tracing. The source is placed
% at target layer, receiver is placed at the free surface or near free
% surface. The receiver position locate above the source position.
%
% Note here the free surface of the velocity model is assume to be 0 (the
% sea-level). Due to the Z axis is vertical down, negative Z values mean
% above the sea-level; and positive values mean below the sea-level (i.e.
% depth). Therefore, the velocity above the sea-level is actually
% undefined. This program currently do not check whether the stations and
% sources are placed above the sea-level (undefined velocity). However, the
% program now just assumes the velocity above the sea-level equals to the
% velocity of the first layer.
%
% Usage:
% input:----------------------------------------------------
% 1 vel: layered velocity model, a vector: nly*1, unit: km/s;
% 2 thk: thickness of each layer, a vector: nly*1, unit: km;
% 3 soup: source position, dimension: nsr*3, each row is a source position,
% column 1: X positon; column 2: Y position; column 3: Z position, unit; km;
% 4 recp: receiver position, dimension: nre*3,each row is a receiver position,
% column 1: X positon; column 2: Y position; column 3: Z position, unit: km.
% output:--------------------------------------------------
% 1 tvt: travel-time, dimension: nsr*nre, row corresponds to source index,
% column corresponds to receiver index; unit: second.
% 2 tkfa: take-off angle, dimension: nsr*nre, row corresponds to source index,
% column corresponds to receiver index; unit: degree; Take-off angle: the
% angle of rays with respect to vertical direction, range: [0 90]; equal to
% incident angle of the rays here;
% 3 rayl: length of the ray, dimension: nsr*nre, row corresponds to source index,
% column corresponds to receiver index, unit: km.

% set tolerate error in comparing distance in km
toler=1e-6;

nsr=size(soup,1); % number of sources
nre=size(recp,1); % number of receivers

% initialize
tvt=zeros(nsr,nre); % travel-times
tkfa=zeros(nsr,nre); % take-off angles
rayl=zeros(nsr,nre); % ray length from source to receiver

% calculate the depth of each layer interface
nly=length(thk); % number of layers
lydp=zeros(nly,1); % do not include free surfaces
lydp(1)=thk(1);
if nly>1
    for ii=2:nly
        lydp(ii)=lydp(ii-1)+thk(ii);
    end
end

% calculate travle-time, take-off angle and ray length for source-reciever pair
for ir=1:nre
    for is=1:nsr
        % calculate the offset of the source-receiver pairs
        srx=recp(ir,1)-soup(is,1); % source-receiver distance in X direction
        sry=recp(ir,2)-soup(is,2); % source-receiver distance in Y direction
        srz=recp(ir,3)-soup(is,3); % source-receiver distance in Z direction
        sroff=sqrt(srx*srx+sry*sry); % source-receiver offset - horizontal diatance (km)
        
        % Calculate which layer the source and receiver belong to
        % Each layer include the bottom interface
        nls=0; nlr=0; % initialization
        for il=1:nly
            if soup(is,3)<=lydp(il)
                nls=il; % source layer
                break
            end
        end
        if nls==0
            error('Wrong source depth at source %d! Exceeding input model depth range!',is);
        end
        for il=1:nly
            if recp(ir,3)<=lydp(il)
                nlr=il; % receiver layer
                break
            end
        end
        if nlr==0
            error('Wrong receiver depth at receiver %d! Exceeding input model depth range!',ir);
        end
        
        % ray-tracing using dichotomy
        if nls~=nlr
            % source and receiver are not at the same layer
            vms=max(vel(1:nls-1)); % find the maximum velocity for the layer above the source
            tkfa1=0; % minimal value of the take-off angle
            if vms<=vel(nls)
                % no total reflection
                tkfa2=90; % maximal value of the take-off angle
            else
                % has higher velocity layer above the source layer
                % take maximum take-off angle which avoid total reflection
                tkfa2=asind(vel(nls)/vms);
            end
            offrd=1; % residual of the calculated source-receiver offset
            szrd=soup(is,3)-lydp(nls-1); % vertical distance in source layer
            fthk=thk;
            fthk(nlr)=lydp(nlr)-recp(ir,3); % vertical distance in reciever layer
            while (offrd>toler)
                tkfat=(tkfa1+tkfa2)/2;
                offtst=szrd*tand(tkfat); % temporary source-receiver offset
                ryltst=sqrt(offtst*offtst+szrd*szrd); % temporary ray-length
                tvttst=ryltst/vel(nls); % temporary travel-time                
                hrp=sind(tkfat)/vel(nls); % horizontal slowness in layered media (s/km)
                for im=nls-1:-1:nlr
                    mcoff=fthk(im)*tand(asind(vel(im)*hrp));
                    offtst=offtst+mcoff;
                    mcryl=sqrt(mcoff*mcoff+fthk(im)*fthk(im));
                    ryltst=ryltst+mcryl;
                    tvttst=tvttst+mcryl/vel(im);
                end
                % reset the bounday of take-off angle according calculation
                if offtst>sroff
                    tkfa2=tkfat; % upper boundary
                else
                    tkfa1=tkfat; % lower boundary
                end
                offrd=abs(sroff-offtst); % calculate the offset residual  
            end
            tvt(is,ir)=tvttst;
            tkfa(is,ir)=tkfat;
            rayl(is,ir)=ryltst;
        else
            % source and receiver are at the same layer
            srdist=sqrt(srx*srx+sry*sry+srz*srz); % source-receiver distance
            tvt(is,ir)=srdist/vel(nlr);
            tkfa(is,ir)=asind(sroff/srdist);
            rayl(is,ir)=srdist;
        end
        
    end
end

end