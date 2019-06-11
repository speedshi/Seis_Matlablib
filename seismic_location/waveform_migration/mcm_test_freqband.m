function [s_pro,n_var]=mcm_test_freqband(fname,stations,mcm,earthquake,search)
% This function is used to show MCM results of different frequency band.
% This can help us determine a suitable frquency band for MCM.
% The MCM results are calculated at specific earthquake location.
%
% INPUT--------------------------------------------------------------------
% fname: file name of the original seismic data;
% stations: matlab structure, contains station information;
% stations.name: cell array, 1*nr, station names;
% stations.north: vector, 1*nr, north coordinates of all stations;
% stations.east: vector, 1*nr, east coordinates of all stations;
% stations.depth: vector, 1*nr, depth coordinates of all stations;
% stations.travelp: P-wave traveltime table, in second, 2D array, ns*nr;
% stations.travels: S-wave traveltime table, in second, 2D array, ns*nr;
% mcm: matlab structure, contains MCM parameters;
% mcm.phasetp: specify seismic phase used for migration, scalar;
% mcm.tpwind: P-phase time window length in second, scalar;
% mcm.tswind: S-phase time window length in second, scalar;
% earthquake: structure, containing the earthquake information;
% earthquake.north: north component of earthquake location, in meter;
% earthquake.east: east component of earthquake location, in meter;
% earthquake.depth: depth component of earthquake location, in meter;
% earthquake.t0: relative origin time of the earthquake, in second, relative
% to the starting time of the seismic data;
% search: matlab structure, describe the imaging area;
% search.soup: source imaging positions, 2D array, ns*3, in meter;
%
% OUTPUT-------------------------------------------------------------------
% s_pro: 2D array, shape: nf*nf, source prominence at different frequency
% band;
% n_var: 2D array, shape: nf*nf, noise variance at different frequency
% band;
%


% read the original seismic data
seismic=read_seis(fname);

% obtain frequency range for testing
fc_s=1; % starting frequency, Hz
fc_e=(seismic.fe)/2-1; % ending frequency, Hz
df=1; % frequency interval, Hz
fc=fc_s:df:fc_e; % the testing frequency points

nf=length(fc); % number of frequency points

% initialize array
s_pro=NaN(nf);
n_var=NaN(nf);

% set the filtering parameters
ffilter.type='bandpass'; % filter type, can be 'low', 'bandpass', 'high', 'stop'
ffilter.order=4; % order of Butterworth filter, for bandpass and bandstop designs are of order 2n

for ifc1=1:nf-1
    for ifc2=ifc1+1:nf
        % set filtering frequency band
        ffilter.freq=[fc(ifc1) fc(ifc2)];
        
        trace=gene_wavetime(seismic,stations,ffilter,[],[],[],[]);
        
        [s_pro(ifc2,ifc1),n_var(ifc2,ifc1)]=mcm_test_para(trace,mcm,search,earthquake,false);
        
    end
end

figure;
fig=imagesc(fc,fc,s_pro);
colormap(jet);colorbar;
set(fig,'AlphaData',~isnan(s_pro)); % set NAN to background color
axis equal tight;
xlabel('f1 (Hz)');ylabel('f2 (Hz)');
title('Source prominence');

figure;
fig=imagesc(fc,fc,n_var);
colormap(flipud(jet));colorbar;
set(fig,'AlphaData',~isnan(n_var)); % set NAN to background color
axis equal tight;
xlabel('f1 (Hz)');ylabel('f2 (Hz)');
title('Noise variance');


end