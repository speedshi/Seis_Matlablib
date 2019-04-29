function model=read_velocity(fname)
% This function is used to read the velocity model. At the present, only
% accept homogeneous model and layered model.
%
% Unit: meter, m/s.
%
% The velocity file format is describe as follow:
% 0                        | float. Reference starting depth (free surface) for the model, can be negtive or positive. 0 means the depth of free surface is the sea-level.
% 100    3500    2000      | float-float-float. Thickness(m)-Vp(m/s)-Vs(m/s) of the layer 1.
% ...                      | float-float-float. Thickness(m)-Vp(m/s)-Vs(m/s) of the layer 2.
%
% NOTE: If there is only one layer, it is a homogeneous model. Negative
% reference starting depth means the free surface is above the sea-level;
% positive reference starting depth means the free surface of the model is
% below the sea-level.
%
% INPUT-----------------------------------------------------
% fname: file name including path of the velocity model;
% OUTPUT----------------------------------------------------
% model: structure, contains velocity model information;
% model.rsd0: reference starting depth of the model, i.e. depth of free surface;
% model.thickness: vector, thickness of each layer;
% model.vp: vector, P-wave velocities of each layer;
% model.vs: vector, S-wave velocities of each layer.

model.rsd0=dlmread(fname,'',[0 0 0 0]); % depth of the free surface of the model

mdps=dlmread(fname,'',1,0); % read in the model parameters

% required parameters
model.thickness=mdps(:,1); % thickness of each layer
model.vp=mdps(:,2); % P-wave velocity of each layer
model.vs=mdps(:,3); % S-wave velocity of each layer

end