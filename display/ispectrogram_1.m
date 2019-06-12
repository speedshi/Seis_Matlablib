function ispectrogram_1(data,dt,name,window)
% This function is used to show the spectrogram of the i-th seismic data;
%
% INPUT--------------------------------------------------------------------
% data: vector, nt*1, the input seismic data;
% dt: scaler, the time sample interval of the data, in second;
% name: string, the name of the station;
% window: scalar, the time window of STFT, in second;


% set default parameter
if nargin < 4
   window=2; % time window of STFT, in second
end

figure;

% show the seismogram
h1=subplot(2,1,1);
plot(data,'k');
title("Seismogram of station: "+name);
axis tight;
set(gca,'XTick',[]);
axsize_1=get(h1,'Position'); % axis position of figure 1

% show the spectrogram
subplot(2,1,2);
show_spectrogram(data,dt,window);
title("Spectrogram of station: "+name);
axsize_2=get(gca, 'Position'); % axis position of figure 2

% make the figure 1 share the same X-axis as figure 2
axsize_1(3)=axsize_2(3); % reset the width of the figure 1
set(h1,'Position',axsize_1);



end