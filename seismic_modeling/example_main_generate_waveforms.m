% This program is used to synthetic waveforms in the
% homogeneous, layered or 3D models.
% The output waveforms are the particle displacements in meter.

% Set input file names of the source, receiver and velocity files
sname="source.dat"; % name of the modeling parameter file for the homogenous and layered modeling
rname="receiver1.dat"; % name of the receiver file for the homogenous and layered modeling
mname_homo="homo.dat"; % file name of the homogeneous velocity model
mname_layer="layer.dat"; % file name of the layered velocity model


% Obtain the source and receiver positions, the source and receiver
% positions are the same for the 3 modelings
[~, ~, ~, source]=rdsourcef(sname); % obtain the source position
stations=rdreceiverf(rname); % obtain the receiver positions
recp = stations.recp;


% For homogeneous model----------------------------------------------------
% show the geometry
plotmodel(sname,rname,mname_homo); % plot the model and monitoring geometries
% generate synthetic data
[ux,uy,uz,dt]=gsynwihl(sname,rname,mname_homo); % synthetic data

% show data
tsmax=1.2; % maximum time in second for showing the data
% maximum of whole dataset is normalized to 1
seisrsdispk(ux,dt); title('ux'); xlim([0 tsmax]);
seisrsdispk(uy,dt); title('uy'); xlim([0 tsmax]);
seisrsdispk(uz,dt); title('uz'); xlim([0 tsmax]);
% maximum of each trace is normailized to 1
seisrsdisp(ux,dt); title('ux'); xlim([0 tsmax]);
seisrsdisp(uy,dt); title('uy'); xlim([0 tsmax]);
seisrsdisp(uz,dt); title('uz'); xlim([0 tsmax]);

% output the data in ESG csv format
outputesgcsv(ux,dt,'ux');
outputesgcsv(uy,dt,'uy');
outputesgcsv(uz,dt,'uz');

% output the data in segy format
wtraces2segy('ux_homo.sgy',ux,dt,source.pos,recp);
wtraces2segy('uy_homo.sgy',uy,dt,source.pos,recp);
wtraces2segy('uz_homo.sgy',uz,dt,source.pos,recp);
%--------------------------------------------------------------------------------------


% For layered model--------------------------------------------------------
% show the geometry
plotmodel(sname,rname,mname_layer); % plot the model and monitoring geometries
% generate synthetic data
[ux,uy,uz,dt]=gsynwihl(sname,rname,mname_layer); % synthetic data

% show data
tsmax=4; % maximum time in second for showing the data
% maximum of whole dataset is normalized to 1
seisrsdispk(ux,dt); title('ux'); xlim([0 tsmax]);
seisrsdispk(uy,dt); title('uy'); xlim([0 tsmax]);
seisrsdispk(uz,dt); title('uz'); xlim([0 tsmax]);
% maximum of each trace is normailized to 1
seisrsdisp(ux,dt); title('ux'); xlim([0 tsmax]);
seisrsdisp(uy,dt); title('uy'); xlim([0 tsmax]);
seisrsdisp(uz,dt); title('uz'); xlim([0 tsmax]);

% save the synthetic data
save ux_layer.mat ux;
save uy_layer.mat uy;
save uz_layer.mat uz;

% output the data in ESG csv format
outputesgcsv(ux,dt,'ux');
outputesgcsv(uy,dt,'uy');
outputesgcsv(uz,dt,'uz');

% output the data in segy format
wtraces2segy('ux_layer.sgy',ux,dt,source.pos,recp);
wtraces2segy('uy_layer.sgy',uy,dt,source.pos,recp);
wtraces2segy('uz_layer.sgy',uz,dt,source.pos,recp);
%--------------------------------------------------------------------------------------


% For 3D modeling----------------------------------------------------------
system('./cdrun.sh','-echo');

% read in the time sample interval for the modeling, in second
fid=fopen('input.dat','r');
xx=textscan(fid,'%f','HeaderLines',9);
fclose(fid);
dt=xx{1};

% read in the generated binary traces
vx=rdfdmodres('vx.rec','input.dat','receiver.dat','single'); % read in the X component
vy=rdfdmodres('vy.rec','input.dat','receiver.dat','single'); % read in the Y component
vz=rdfdmodres('vz.rec','input.dat','receiver.dat','single'); % read in the Z component

% show data
tsmax=4; % maximum time in second for showing the data
% maximum of whole dataset is normalized to 1
seisrsdispk(vx,dt); title('vx'); xlim([0 tsmax]);
seisrsdispk(vy,dt); title('vy'); xlim([0 tsmax]);
seisrsdispk(vz,dt); title('vz'); xlim([0 tsmax]);
% maximum of each trace is normailized to 1
seisrsdisp(vx,dt); title('vx'); xlim([0 tsmax]);
seisrsdisp(vy,dt); title('vy'); xlim([0 tsmax]);
seisrsdisp(vz,dt); title('vz'); xlim([0 tsmax]);

% transform from velocity to displacement
ux=trvel2dis(vx,dt);
uy=trvel2dis(vy,dt);
uz=trvel2dis(vz,dt);

% output to segy file
wtraces2segy('ux_3D.sgy',ux,dt,source.pos,recp);
wtraces2segy('uy_3D.sgy',uy,dt,source.pos,recp);
wtraces2segy('uz_3D.sgy',uz,dt,source.pos,recp);
%---------------------------------------------------------------------------------------
