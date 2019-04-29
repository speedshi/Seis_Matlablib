function lydp=callyifdp(thickness,rsd0)
% This function is used to calculate the layer depth of each interface for
% the layered model (including the free surface).
% Unit: meter.
% INPUT--------------------------------------------------------------------
% rsd0: the depth of the free surface (top of the first layer),
% 0-->sea-level, (+)positive-->below sea-level, (-)negative-->above sea-level;
% thickness: the thickness of each layer, a vector.
% OUTPUT-------------------------------------------------------------------
% lydp: the depths of each layer interface (lay boundaries) including the free surface.

% set default value
if nargin < 2
    rsd0=0; % default free surface of the model is at the sea-level
end

nly=max(size(thickness))+1; % the number of interfaces (layer boundaries)
lydp=zeros(nly,1);

% depth of each layer surface, including the free surface
lydp(1)=rsd0;
for ii=2:nly
    lydp(ii)=lydp(ii-1)+thickness(ii-1);
end

end