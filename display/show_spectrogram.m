function show_spectrogram(data,dt,window)
% This function is used to show the spectrogram of the input seismic data;
%
% INPUT--------------------------------------------------------------------
% data: seismic data, vector, shape:nt*1;
% dt: time sampling interval of the input seismic data, in second;
% window: the time window of the Short-time Fourier transform, in second;
%


fs=1/dt; % sampling frequency in Hz

window=window/dt; % window of STFT in points

spectrogram(data,window,[],[],fs,'yaxis');
colormap(jet);
%ax=gca; ax.YScale='log';


end