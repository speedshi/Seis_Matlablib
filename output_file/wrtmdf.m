function wrtmdf(thk,vp,vs,den,qp,qs,fname)
% This function is used to write the elastic parameter file of the model.
% The file is used for 'fk' model input.
% Usage: 
% input: thk----thickness of each layer (m);
%           vp----P-wave velocity of each layer (m/s);
%           vs----S-wave velocity of each layer (m/s);
%           den---density of each layer (kg/m^3);
%           qp----absorbing value of P-wave of each layer;
%           qs----absorbing value of S-wave of each layer;
%           fname--the name of the output file.
% 'thk', 'vp', 'vs', 'den', 'qp' and 'qs' must be a vector and have the same length.
% output: a file whose name is 'model', (saved at '../data/').
% For 'fk' input, the file format of the 'model' file is:
% thickness | vs | vp | density | Qs | Qp
% the unit of thickness is km; the unit of velocity is km/s; the unit of
% density is g/cm^3.

if nargin<7
    fname='model'; % set default name
end

% transfer unit: m -> km, m/s -> km/s, kg/m^3 -> g/cm^3
tthk=thk/1000; % m -> km
tvp=vp/1000; % m/s -> km/s
tvs=vs/1000; % m/s -> km/s
tden=den/1000; % kg/m^3 -> g/cm^3

fnm=fname; % save data file into a specified folder
fid1=fopen(fnm,'wt');
nl=max(size(thk)); % number of layers

for in=1:nl
    fprintf(fid1,'%.4f  %.4f  %.4f  %.4f  %.2f  %.2f  \n',tthk(in),tvs(in),tvp(in),tden(in),qs(in),qp(in));
end
fclose(fid1);
end