function data=cal_charfunc(data,tp)
% This function is used to calculate the characteristic function for
% migration.
%
% INPUT--------------------------------------------------------------------
% data: the input seismic data, 2D array, shape: nt*nre;
% tp: specify the type of the characteristic function, scalar;
%
% OUTPUT:------------------------------------------------------------------
% data: the calculated characteristic function, shape: nt*nre;

switch tp
    case 0
        % original data
        fprintf('Use original waveform data as the characteristic function.\n');
    case 1
        % envelope
        fprintf('Use envelope as the characteristic function.\n');
        data=abs(hilbert(data));
    case 2
        % absolute values
        fprintf('Use absolute value as the characteristic function.\n');
        data=abs(data);
    case 3
        % non-negative values
        fprintf('Use non-negative value as the characteristic function.\n');
        data(data<0)=0;
    case 4
        % square values
        fprintf('Use square value as the characteristic function.\n');
        data=data.^2;
    otherwise
        error('Incorrect input for mcm.cfuntp!');
end





end