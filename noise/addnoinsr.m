function data=addnoinsr(signal,noise,nsr)
% This function is used to add noise into the signal according to specified
% noise-to-signal ratio.
% Input:----------------------------------------------------
% signal: signal data;
% noise: noise data which has the same dimension with signal;
% nsr: noise-to-signal ratio, determined by the ratio of the maximum
% amplitude of the noise data and signal data.
% Output:--------------------------------------------------
% data: data after adding noise.

r=nsr*max(abs(signal(:)))/max(abs(noise(:))); % calculate the ratio to apply on the noise data
data=signal+r*noise; % adjust the level of noise data and add noise into the signal data

end