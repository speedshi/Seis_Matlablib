function [vp_rms,vs_rms]=equihomo(thk,vp,vs,den,qp,qs,fname)
% This function is used to calculate the equivalent root-mean-square
% velocity of the layered media. Then output the equivalent homogeneous
% model for 'fk' model input.
% Usage: 
% input: thk----thickness of each layer (m);
%           vp----P-wave velocity of each layer (m/s);
%           vs----S-wave velocity of each layer (m/s);
%           den---density of each layer (kg/m^3);
%           qp----absorbing value of P-wave of each layer;
%           qs----absorbing value of S-wave of each layer;
%           fname--the name of the output file.
% 'thk', 'vp', 'vs', 'den', 'qp' and 'qs' must be a vector and have the same length.
% OUTPUT (a file whose default name is 'modelhm'.)
% vp_rms: root-mean-square of P-wave for the input laryered model;
% vs_rms: root-mean-square of S-wave for the input laryered model.
% For 'fk' input, the file format of the 'model' file is:
% thickness | vs | vp | density | Qs | Qp
% the unit of thickness is km; the unit of velocity is km/s; the unit of
% density is g/cm^3.

if nargin<7
    fname='modelhm'; % set default name
end

% transfer unit: m -> km, m/s -> km/s, kg/m^3 -> g/cm^3
tthk=thk/1000; % m -> km
tvp=vp/1000; % m/s -> km/s
tvs=vs/1000; % m/s -> km/s
tden=den/1000; % kg/m^3 -> g/cm^3

nl=max(size(thk)); % number of layers
if nl>1
    tvp2sm=0;
    tvs2sm=0;
    tpism=0;
    tsism=0;
    for il=1:nl-1
        % use the top nly-1 layers to calculate the V_rms, the last layer
        % is half space.
        tpi=tthk(il)/tvp(il); % vertical P-wave travel-time for each layer
        tsi=tthk(il)/tvs(il); % vertical S-wave travel-time for each layer
        tvp2sm=tvp2sm+tpi*tvp(il)*tvp(il); % sum(t_i*vp_i^2)
        tvs2sm=tvs2sm+tsi*tvs(il)*tvs(il); % sum(t_i*vs_i^2)
        tpism=tpism+tpi; % sum(tp_i)
        tsism=tsism+tsi; % sum(ts_i)
    end
    vp_rms=sqrt(tvp2sm/tpism); % root-mean-square of P-wave
    vs_rms=sqrt(tvs2sm/tsism); % root-mean-square of S-wave
    den_ave=mean(tden(1:nl-1)); % average value of density
    qp_ave=mean(qp(1:nl-1)); % average value of Qp
    qs_ave=mean(qs(1:nl-1)); % average value of Qs    
else
    % only one layer, no need to calculate
    vp_rms=tvp;
    vs_rms=tvs;
    den_ave=tden;
    qp_ave=qp;
    qs_ave=qs;    
end

% write to file
fnm=fname; % save data file into a specified folder
fid1=fopen(fnm,'wt');
fprintf(fid1,'%.4f  %.4f  %.4f  %.4f  %.2f  %.2f  \n',0,vs_rms,vp_rms,den_ave,qs_ave,qp_ave);
fclose(fid1);
end