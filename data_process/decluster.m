function [xx_new, yy_new] = decluster(xx,yy,xlim)
% This function is used to decluster the input dataset.
% Timing or indexing ('xx') of the input signals/events within the set
% threshold 'xlim' are viewed as the same signals/events, thus only keep
% the one with the maximum Y values ('yy').
%
% INPUT:
%       xx: the timing or indexing of the input signals/events;
%       yy: the feature values (such as CC) of the input signals/events;
%       xlim: the timing or indexing threshold, smaller then (<=) this value
%             will be viewed as the same signals/events;
%
% OUTPUT:
%       xx_new: the timing or indexing of the declustered signals/events;
%       yy_new: the feature values of the declustered signals/events;
%


NN = length(xx);  % total number of input events/signals

% initial the outputs
xx_new = [];
yy_new = [];

xx_temp = xx(1);
yy_temp = yy(1);

for ii = 2:NN
    if (xx(ii) - xx(ii-1)) <= xlim
        % belong to the same signals/events
        if yy(ii) > yy_temp
            xx_temp = xx(ii);
            yy_temp = yy(ii);
        end
        
    else
        % are different signals/events
        xx_new(end+1) = xx_temp;
        yy_new(end+1) = yy_temp;
        
        xx_temp = xx(ii);
        yy_temp = yy(ii);
    end
end

% For the last one, need to add into the final results; since when ii=NN, 
% the last signal/event are not added into the outputs yet.
xx_new(end+1) = xx_temp;
yy_new(end+1) = yy_temp;

% convert to vector
xx_new = xx_new(:);
yy_new = yy_new(:);

end