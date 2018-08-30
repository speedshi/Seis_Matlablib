function datans=pnoise(data,snr,seed)
% This function is used to add Gaussian random noise into the data.
% Input:***********************************
% (1) data: the input noise free data;
% (2) snr: signal-noise-ratio in terms of energy ratio;
% (3) seed: seeds the random number generator.
% Output:**********************************
% (1) datans: data after adding noise.

if nargin==3
    rng(seed); % set the seed for randn
end

s_ener= norm(data(:))^2; % the energy of original data
zao=randn(size(data)); % generate white Gaussian noise
zao_ener=norm(zao(:))^2; % the initial energy of the noise
factor=sqrt(s_ener/zao_ener/snr); % calculate the scale factor for noise
noise=factor*zao; % apply the scale factor, keep 'SNR' consistent
datans=data+noise; % add noise on the data