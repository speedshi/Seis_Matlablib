function ispectrogram(seismic,i)
% This function is used to show the spectrogram of the i-th seismic data;
%
% INPUT--------------------------------------------------------------------
% seismic: structure, contains seismic data and metadata of each station;
% seismic.name: cell array, 1*ns, contains the name of each station;
% seismic.fe: scaler, the sampling frequency of the data, in Hz;
% seismic.dt: scaler, the time sample interval of the data, in second;
% seismic.t0: matlab datetime, the origin time of the seismic data;
% seismic.data: 2D array, ns*nt, contains seismic data;
% seismic.network: string, the name of the network;
% seismic.component: character, the name of the data component (usually N, E or Z);
% i: scalar, specifiy to show the spectrogrma of which station;


% set default parameter
if nargin <2
   i=1;
end

figure;

% show the seismogram
h1=subplot(2,1,1);
plot(seismic.data(i,:),'k');
title("Seismogram of station: "+seismic.name{i});
axis tight;
set(gca,'XTick',[]);
axsize_1=get(h1,'Position'); % axis position of figure 1

% show the spectrogram
window=2; % time window of STFT, in second
subplot(2,1,2);
show_spectrogram(seismic.data(i,:),seismic.dt,window);
title("Spectrogram of station: "+seismic.name{i});
axsize_2=get(gca, 'Position'); % axis position of figure 2

% make the figure 1 share the same X-axis as figure 2
axsize_1(3)=axsize_2(3); % reset the width of the figure 1
set(h1,'Position',axsize_1);



end