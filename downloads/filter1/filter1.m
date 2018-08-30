function [yfilt,filtb,filta] = filter1(filtertype,y,varargin) 
% filter1 performs frequency or wavelength filtering on a 1D array 
% using zero-phase Butterworth filtering. 
% 
%% Syntax
% 
%  yfilt = filter1(filtertype,y,'fc',Fc)
%  yfilt = filter1(filtertype,y,'lambdac',lambdac)
%  yfilt = filter1(...,'fs',Fs)
%  yfilt = filter1(...,'x',x)
%  yfilt = filter1(...,'Ts',Ts)
%  yfilt = filter1(...,'order',FilterOrder)
%  [yfilt,filtb,filta] = filter1(...)
% 
%% Description
% 
% yfilt = filter1(filtertype,y,'fc',Fc) filters 1D signal y
% using a specified filtertype and cutoff frequency Fc. For 
% high-pass or low-pass filters Fc must be a scalar. For band-
% pass and band-stop filters Fc must be a two-element array. The 
% filtertype can be 
% 
%   * 'hp' high-pass with scalar cutoff frequency Fc
%   * 'lp' low-pass with scalar cutoff frequency Fc
%   * 'bp' band-pass with two-element cutoff frequencies Fc
%   * 'bs' band-stop with two-element cutoff frequencies Fc 
% 
% yfilt = filter1(filtertype,y,'lambdac',lambdac) specifies cutoff
% wavelength(s) rather than cutoff frequencies.  This syntax assumes
% lambda = 1/f. 
% 
% yfilt = filter1(...,'fs',Fs) specifies a sampling frequency Fs. 
% If neither 'fs', 'x', nor 'Ts' are specified, Fs = 1 is assumed.    
% 
% yfilt = filter1(...,'x',x) specifies a vector of  monotonically-
% increasing, equally-spaced sampling times or x locations corresponding
% to y, which is used to determine sampling frequency. If neither 'fs', 
% 'x', nor 'Ts' are specified, Fs = 1 is assumed.  
% 
% yfilt = filter1(...,'Ts',Ts) specifies a sampling period or sampling distance
% such that Fs = 1/Ts. If neither 'fs', 'x', nor 'Ts' are specified, 
% Fs = 1 is assumed.    
% 
% yfilt = filter1(...,'order',FilterOrder) specifies the order (sometimes 
% called rolloff) of the Butterworth filter. If unspecified, FilterOrder = 1 is assumed. 
% 
% [yfilt,filtb,filta] = filter1(...) also returns the filter numerator 
% filta and denominator filtb. 
% 
%% Example 
% % For this example we'll use the built-in train whistle example file and we'll add 
% % a little gaussian random noise to make things interesting.  
% 
%     load train 
%     y = y+0.1*randn(size(y)); 
% 
% % High-pass filter the train whistle, keeping only frequencies above 750 Hz:  
% 
%     yhp = filter1('hp',y,'fs',Fs,'fc',750); 
% 
% % Lowpass filter the train whistle, keeping only frequencies below 1100 Hz
% % and specify a sharp 5th order rolloff. Also reference y to a time vector
% % instead of the sampling rate we specified above: 
% 
%     t = (0:length(y)-1)/Fs; 
%     ylp = filter1('lp',y,'x',t,'fc',1100,'order',5); 
%     
% % Bandstop filter the train whistle from 750 Hz to 1000 Hz to eliminate the 886 Hz 
% % middle frequency.      
%   
%     ybs = filter1('bs',y,'fs',Fs,'fc',[750 1000],'order',5); 
% 
% % Use plotpsd (available on File Exchange) to show original and
% % bandstop-filtered signals: 
% 
%     plotpsd(y,Fs,'b')
%     hold on
%     plotpsd(ybs,Fs,'r')
% 
%% Author Info
% The filter1 function was written by Chad A. Greene of the University of
% Texas at Austin's Institute for Geophysics (UTIG), October 2015. 
% http://www.chadagreene.com
% 
% See also butter and filtfilt. 

%% Initial error checks: 

assert(license('test','signal_toolbox')==1,'The filter1 function requires Matlab''s Signal Processing Toolbox.')
assert(nargin>3,'Not enough input arguments.') 
assert(sum(strcmpi({'hp';'lp';'bp';'bs';'high';'low';'bandpass';'stop'},filtertype))==1,'Filter type must be ''hp'', ''lp'', or ''bp''.'),
assert(sum([strcmpi(varargin,'fc') strcmpi(varargin,'lambdac')])==1,'Must declare a cutoff frequency (or frequencies) ''fc'', or cutoff wavelength(s) ''lambdac''.')
assert(isvector(y)==1,'Input y must be a vector.') 

%% Define defaults: 

order = 1; 
Fs = 1; 

%% Parse Inputs: 

% Replace filtertype string if necessary: 
filtertype = strrep(filtertype,'hp','high'); 
filtertype = strrep(filtertype,'lp','low'); 
filtertype = strrep(filtertype,'bp','bandpass'); 
filtertype = strrep(filtertype,'bs','stop'); 

% Is a sampling frequency defined? 
tmp = strcmpi(varargin,'fs'); 
if any(tmp) 
    Fs = varargin{find(tmp)+1}; 
end

% Define sampling period: 
tmp2 = strcmp(varargin,'Ts'); 
if any(tmp2)
    Fs = 1/varargin{find(tmp2)+1}; 
end

% Define sampling vector: 
tmp3 = strcmp(varargin,'x'); 
if any(tmp3) 
    x = varargin{find(tmp3)+1}; 
    assert(isvector(x)==1,'Input x must be a vector.')
    assert(length(x)==length(y),'Dimensions of input vector x must match dimensions of input signal vector y.') 
    
    Ts = unique(diff(x));
    assert(all([isscalar(Ts) isfinite(Ts)])==1,'Input vector x must be equally spaced.')
    Fs = 1/Ts; 
end

% Make sure user didn't try to define a sampling frequency AND a sampling period: 
assert(any(tmp)+any(tmp2)+any(tmp3)<2,'I am confused. It looks like you have attempted to define a sampling frequency and a sampling period.  Check inputs of filter1.') 


% Cutoff Frequency: 
tmp = strcmpi(varargin,'fc'); 
if any(tmp) 
    cutoff_freqs = varargin{find(tmp)+1}; 
end

tmp2 = strncmpi(varargin,'lambdac',4); 
if any(tmp2) 
    cutoff_freqs = 1./varargin{find(tmp2)+1}; 
end
assert(any(tmp)+any(tmp2)<2,'I am confused. It looks like you have attempted to define a cutoff frequency and a cutoff period.  Check inputs of filter1.') 

% Filter order: 
tmp = strncmpi(varargin,'order',3); 
if any(tmp) 
    order = varargin{find(tmp)+1}; 
end

%% Error checks on inputs: 

assert(isscalar(Fs)==1,'Input error: Undefined sampling frequency or period.')

switch filtertype 
    case {'low','high'} 
        assert(isscalar(cutoff_freqs)==1,'Low-pass and High-pass filters require a scalar cutoff frequency.') 
        
    case {'stop','bandpass'} 
        assert(numel(cutoff_freqs)==2,'Bandpass and bandstop filters require a low and high frequency.') 
        cutoff_freqs = sort(cutoff_freqs); 
        
    otherwise
        error('Unrecognized filter type.') 
end
      
%% Construct filter: 

nyquist_freq = Fs/2;             % Nyquist frequency
Wn=cutoff_freqs/nyquist_freq;    % non-dimensional frequency
[filtb,filta]=butter(order,Wn,filtertype); % construct the filter
yfilt=filtfilt(filtb,filta,y);   % filter the data with zero phase 

end