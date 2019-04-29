function data=fdmodres(fname,nrec,nt,prec)
% This function is used to read the output traces of the fdmodeling
% software. Note the output traces is in binary format, and the storage
% sequence is 'nrec-nt'. However, the final output data is in nt*nrec.
% INPUT----------------------------------------------------
% fname: name of the trace file including the path;
% nrec: total number of traces;
% nt: total number of time samples;
% prec: data precision of the traces, 'single' or 'double'.
% OUTPUT-------------------------------------------------
% data: read in data, nt*nrec.

if nargin<4
    prec='single'; % default is in single precision
end

fid=fopen(fname,'r');
data=fread(fid,[nrec nt],prec);
fclose(fid);

data=data'; % note the sequence is nt*nrec

end