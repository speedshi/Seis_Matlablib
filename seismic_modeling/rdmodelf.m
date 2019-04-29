function [nl,rsd0,model]=rdmodelf(fname)
% Read the model parameter file.
% INPUT:-------------------------------------------------
% fname: file name of the model parameter file.
% OUTPUT:----------------------------------------------
% nl: number of layers of the model;
% rsd0: reference starting depth of the model, i.e. depth of free surface;
% model: structure array containing model parameters;
% model.thickness: thickness of each layer, nl*1;
% model.vp: P-wave velocity of each layer, nl*1;
% model.vs: S-wave velocity of each layer, nl*1;
% model.den: density of each layer, nl*1;
% model.qp: quality factor (attenuation) of P-wave, nl*1, optional;
% model.qs: quality factor (attenuation) of S-wave, nl*1, optional.


if nargin<1
    fname='model.dat';
end

rsd0=dlmread(fname,'',[0 0 0 0]); % depth of the free surface of the model

mdps=dlmread(fname,'',1,0); % read in the model parameters

[nl,npr]=size(mdps); % number of layers and elastic parameters of the model

% required parameters
model.thickness=mdps(:,1); % thickness of each layer
model.vp=mdps(:,2); % P-wave velocity of each layer
model.vs=mdps(:,3); % S-wave velocity of each layer
model.den=mdps(:,4); % density of each layer

% optional parameters
if npr==4
    model.qp=[];
    model.qs=[];
elseif npr==5
    model.qp=mdps(:,5);
    model.qs=[];
elseif npr==6
    model.qp=mdps(:,5); % quality factor of P-wave
    model.qs=mdps(:,6); % quality factor of S-wave
else
    error('Input model parameter file does not have a correct format.');
end

end