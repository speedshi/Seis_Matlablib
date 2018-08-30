function datans=pnoisem(data,ncp,seed)
% This function is used to add Gaussian random noise into the data.
% The noise is expressed as the ratio of maximum amplitude between noise
% and signal. ncp=noise_max/signal_max.
% Input:***********************************
% (1) data: the input signal data;
% (2) ncp: noise level relative to the maximum amplitude;
% (3) seed: seeds the random number generator.
% Output:**********************************
% (1) datans: data after adding noise.

if nargin==3
    rng(seed); % set the seed for randn
end

zao=randn(size(data)); % generate white Gaussian random noise
dmv=max(abs(data(:))); % maximum amplitude of signal
zmv=max(abs(zao(:))); % maixmum amplitude of noise
factor=ncp*dmv/zmv; % calculate scale factor for noise according to ncp
noise=factor*zao; % apply the scale factor
datans=data+noise; % add noise on the data