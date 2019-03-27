function idata=intpmodel3(data,lx,ly,lz,dx,dy,dz)
% This function is used to interpolate the input data according to the
% estimated model length and grid interval in each dimension.
% Need to make sure lx/dx, ly/dy and lz/dz are all integers.
% NOTE the dimension of the input 3D 'data' is X-Y-Z; first dimension is X,
% second dimension is Y, third dimension is Z.
% INPUT---------------------------------------------------------------
% data: 3D original data, nx*ny*nz, note it is in X-Y-Z dimension order;
% lx: estimated length in the X direction, in meter;
% ly: estimated length in the Y direction, in meter;
% lz: estimated length in the Z direction, in meter;
% dx: grid interval in the X direction, in meter;
% dy: grid interval in the Y direction, in meter;
% dz: grid interval in the Z direction, in meter.
% OUTPUT-----------------------------------------------------------
% idata: output data after interpolation.

[nx,ny,nz]=size(data); % the original size of the data

enx=lx/dx+1; % estimated grid size in the X direction
eny=ly/dy+1; % estimated grid size in the Y direction
enz=lz/dz+1; % estimated grid size in the Z direction

% Interpolation
xx=linspace(1,enx,nx);
yy=linspace(1,eny,ny);
zz=linspace(1,enz,nz);
[X,Y,Z] = meshgrid(xx,yy,zz);
[Xq,Yq,Zq] = meshgrid(1:enx,1:eny,1:enz);
rdata=permute(data,[2 1 3]); % rearrange dimensions of the input array to "Y-X-Z"
idata=interp3(X,Y,Z,rdata,Xq,Yq,Zq);
idata=permute(idata,[2 1 3]); % rearrange back to "X-Y-Z"

end