function [nr,recp]=rdreceiverf(fname)
% Read the receiver file.
% INPUT:-------------------------------------------------
% fname: file name of the receiver file.
%OUTPUT:-----------------------------------------------
% nr: number of recerivers
% recp: X-Y-Z positions of receivers, nr*3.

if nargin<1
    fname='receiver.dat';
end

rloc=importdata(fname,' ',1);
nr=size(rloc.data,1); % number of receivers
recp=rloc.data; % X-Y-Z locations

end