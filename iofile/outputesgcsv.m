function outputesgcsv(data,dt,fname,t0)
% This function is used to output the data in ESG csv data format.
% INPUT:-------------------------------------------------
% data: the seismic waveform data,nt*nr;
% dt: time sampling rate, in second;
% t0: starting time of the record, in second, time are calculated using: t0+(n-1)*dt;
% fname: output file name.

if nargin<4
    t0=0; % set default starting time as 0 second.
end

[nt,nr]=size(data);

fid=fopen(fname,'wt');

ctime=char(datetime('now'));

fprintf(fid,'%s%s%sE  ESG   %s\n',ctime(4:6),ctime(1:2),ctime(13:20),ctime(8:11));
fprintf(fid,'PAD #:       %s\n',fname);
fprintf(fid,' NUMBER OF SAMPLES PER CHANNEL=%8d\n',nt);
fprintf(fid,' NUMBER OF TRIAXIAL CHANNELS..=%8d\n',nr);
fprintf(fid,' NUMBER OF UNIAXIAL CHANNELS..=%8d\n',0);
fprintf(fid,' NUMBER OF CHANNELS EXTRACTED.=%8d\n',nr);
fprintf(fid,' SAMPLING RATE................=%8d\n',1/dt);
fprintf(fid,'   MSEC,'); % note the output time unit is millisecond

for ir=1:nr
    fprintf(fid,'     CH%3d,    ',ir);
end
fprintf(fid,'\n');

for it=1:nt
    fprintf(fid,'%7.2f,',(t0+(it-1)*dt)*1000);
    for ir=1:nr
        fprintf(fid,'%14.5e,',data(it,ir));
    end
    fprintf(fid,'\n');
end

fclose(fid);


end