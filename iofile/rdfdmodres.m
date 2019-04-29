function data=rdfdmodres(fname,inputname,recfname,prec)
% This function is used to read the output traces of the fdmodeling
% software. Note the output traces is in binary format, and the storage
% sequence is 'nrec-nt'. However, the final output data is in nt*nrec.
% Need to keep consistant with the format of the input file 'input.dat' and
% the receiver file 'receiver.dat' of the fdmodeling software.
% INPUT----------------------------------------------------
% fname: name of the trace file including the path;
% inputname: name of the input file of the fdmodeling software;
% recfname: name of the receiver file of the fdmodeling software;
% prec: data precision of the traces, 'single' or 'double'.
% OUTPUT-------------------------------------------------
% data: read in data, nt*nrec.

if nargin<4
    prec='single'; % default is in single precision
end

% read in the total number of time samples for the modeling
fid=fopen(inputname,'r');
xx=textscan(fid,'%d','HeaderLines',8);
fclose(fid);
nt=double(xx{1});

% read in the total numver of receivers for the modeling
rloc=importdata(recfname,' ',1);
nrec=size(rloc.data,1); % number of receivers


fid=fopen(fname,'r');
data=fread(fid,[nrec nt],prec);
fclose(fid);

data=data'; % note the sequence is nt*nrec

end