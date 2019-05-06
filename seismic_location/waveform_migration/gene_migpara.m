function gene_migpara(mcm,fname)
% This function is used to generate the text file for MCM parameters.
%
% INPUT--------------------------------------------------------------------
% mcm: matlab structure, contains all the required MCM parameters;
% fname: file name of the output text file.
% OUTPUT-------------------------------------------------------------------
% migpara.dat: text file for MCM parameters.

folder='data'; % name of the folder where output data are stored

% check if the output folder exists, if not, then create it
if ~exist(folder,'dir')
    mkdir(folder);
end

% set default value
if nargin==1
   fname='migpara.dat'; 
end

fname=['./' folder '/' fname]; % including the folder

fid=fopen(fname,'wt');

fprintf(fid,'%d\n',mcm.migtp);
fprintf(fid,'%d\n',mcm.phasetp);
fprintf(fid,'%d\n',mcm.cfuntp);
fprintf(fid,'%d\n',mcm.nre);
fprintf(fid,'%d\n',mcm.nsr);
fprintf(fid,'%s\n',mcm.dfname);
fprintf(fid,'%f\n',mcm.dt);
fprintf(fid,'%f\n',mcm.tdatal);
fprintf(fid,'%f\n',mcm.tpwind);
fprintf(fid,'%f\n',mcm.tswind);
fprintf(fid,'%f\n',mcm.dt0);
fprintf(fid,'%f\n',mcm.vthrd);
fprintf(fid,'%d\n',mcm.mcmdim);
fprintf(fid,'%f\n',mcm.spaclim);
fprintf(fid,'%f\n',mcm.timelim);
fprintf(fid,'%d\n',mcm.nssot);

fclose(fid);

end