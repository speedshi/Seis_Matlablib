function [wzc,wnc,wec,fzc,fnc,fec]=wavefk(soup,recp,recz,smt,sm0,st0,dt,ENt)
% This function is used to call 'fk' to do full waveform modeling in layered
% media. The convention of coordinate systerm follows the definition in Aki &
% Richards 2002, Fig 4.20 and BOX 4.4. Right-hand coordinate systerm in
% Cartesian coordinate, with X-North, Y-East, Z-Vertical down.
% The unit of length in this program are all in 'km'.
% For seismic moment tensor, the order of each component are [mxx  mxy  mxz  myy  myz  mzz].
% The azimuth is measured based on the source point clockwise from North (X axis) to geophone point.
% All the input and output data are stored and saved at folder '../data'.
% INPUT-------------------------------------------------------------------------------
% soup: source position. Each row represents a three-dimensional source
% position (column 1: X, column 2: Y, column 3: Z), nsr*3;
% recp: surface geophone positioin. Each row represents a two-dimensional
% geophone position (column 1: X, column 2: Y),nre*2, all geophones share the same depth;
% recz: depth of the goephones, scaler;
% smt: moment tensors of each source, nsr*6, order: [mxx  mxy  mxz  myy myz  mzz];
% sm0: seismic moment of earch source, nsr*1, unit: dyne-cm;
% st0: origin times (onset time) of each source relative to the
% starting time of the recorded data, in second, nsr*1, should be positive;
% dt: the time interval of the recorded data, in second, scaler;
% ENt: number of time samples for the modeling data, scaler;
% OUTPUT-------------------------------------------------------------------------------
% (2) wzc: Z component full waveform data of different sources, ENt*nre*nsr;
% (3) wnc: N component full waveform data of different sources, ENt*nre*nsr;
% (4) wec: E component full waveform data of different sources, ENt*nre*nsr;
% (5) fzc: Z component full waveform data with all the sources added together, ENt*nre;
% (6) fnc: N component full waveform data with all the sources added together, ENt*nre;
% (7) fec: E component full waveform data with all the sources added together, ENt*nre.
%*******************************************************************************
% Required file:
% (2) 'model': text file used as the 'fk' input file;
% (3) 'fk.pl': perl script file used to call 'fk' to do waveform modelling;
% (4) 'fk': executive file of 'fk', for modelling using reflectivity method;
% (5) 'trav': executive file for calculating arrival times;
% (6) 'sachd': executive file for adding and modifying sac data head;
% (7) 'st_fk': executive file of 'fk', for calculating static displacement;
% (8) 'syn': executive file for synthesizing full waveform data using calculated Green's function;
% (9) 'readsacp': executive file for transferring sac data into binary data without data head.
%*******************************************************************************

% parameters settings--------------------------------------------------------------
dpre=5; % the precision of distance in km, expressed as the number after decimal point, keep consistent with 'fk'
%-----------------------------------------------------------------------------------------------------------

pp=nextpow2(ENt);
Nt=2^pp; % the total time points of the recorded data for 'fk', should be 2^N.

nsr=size(soup,1); % number of sources
nre=size(recp,1); % number of surface geophones
stf=zeros(nsr,Nt); % source time function for each source
% calculate source time function, here we assume the Ricker wavelet apply
% for all the sources, and frequency is 40 Hz. Note the time delay for
% Ricker wavelet is set as 1.1/freq by default.
freq=40;
for is=1:nsr
    stf(is,:)=rickerw(freq,dt,Nt);
end

% calculate the offset and azimuth of the source-receiver pairs
sroff=zeros(nsr,nre); % offsets (km)
srazi=zeros(nsr,nre); % azimuth (degree)
for is=1:nsr
    for ir=1:nre
        % calculate the source-receiver offset
        srx=recp(ir,1)-soup(is,1); sry=recp(ir,2)-soup(is,2);
        sroff(is,ir)=sqrt(srx^2+sry^2);
        % calculate the azimuth angle of the geophone (degree)
        if sroff(is,ir)<=eps
            srazi(is,ir)=0; % geophone is just placed above the source point
        else
            if sry>=0
                srazi(is,ir)=acosd(srx/sroff(is,ir));
            else
                srazi(is,ir)=360-acosd(srx/sroff(is,ir));
            end
        end
    end
end

% set precision of source depth and source-receiver offsets
soup(:,3)=round(soup(:,3),dpre); % set precision of the source depth, need to set precision before use 'unique' to pick and order the depth
sroff=round(sroff,dpre); % set precision of the source-receiver offsets, need to set precision before use 'unique' to pick and order the offsets

% use 'fk' to generate Green's function
sdp=unique(soup(:,3)); % find out different depth of sources, and also sort the depth in ascending order
nsdp=max(size(sdp)); % number of different source depth
status1=zeros(nsdp,1);status2=zeros(nsdp,1); % used to save the 'fk' running status
cmdout1=cell(nsdp,1);cmdout2=cell(nsdp,1); % used to save the output of 'fk'
sftp=sprintf('%%.%df',dpre); % get string for setting precision in 'num2str'
sftps=sprintf('%%.%df  ',dpre); % get string for setting precision in 'num2str', note the space, space is important to separate the numbers
sdt=num2str(dt); % transfer time interval to string, note here we adopt the default precision, should be enough
srez=num2str(recz,sftp); % transfer geophone depth to string
for ii=1:nsdp
    sid=soup(:,3)==sdp(ii); % index of this particular source depth, logical array
    srf1=sroff(sid,:); % pick out the source-receiver offset of this depth, merge them
    srf1=srf1(:)';
    srf=unique(srf1); % only save the different source-reveiver offsets, delete the same offsets of this depth
    ssdp=num2str(sdp(ii),sftp); % transfer source depth to string
    ssrf=num2str(srf,sftps); % transfer source-receiver depth to string, note the 'blank space' is very important for separating the different numbers
    sfk1=sprintf('fk.pl -Mmodel/%s  -H0/0  -N%d/%s/1/0.1/0.3  -P0/1  -R%s  -S0  %s',ssdp,Nt,sdt,srez,ssrf); % for explosive source
    sfk2=sprintf('fk.pl -Mmodel/%s  -H0/0  -N%d/%s/1/0.1/0.3  -P0/1  -R%s  -S2  %s',ssdp,Nt,sdt,srez,ssrf); % for double-couple source
    [status1(ii), cmdout1{ii}]=system(sfk1,'-echo'); % call 'fk' to generate explosive green's function
    [status2(ii), cmdout2{ii}]=system(sfk2,'-echo'); % call 'fk' to generate double-couple green's function
end

% use the generated Green's function to synthesize full waveform data
system('mkdir sdata'); % create a folder to save the synthetic data
for is=1:nsr
    for ir=1:nre
        mxx=num2str(smt(is,1));  % transfer the seismic moment mxx to string, use default precision (4 number after decimal point), should be enough
        mxy=num2str(smt(is,2));  % transfer the seismic moment mxy to string, use default precision (4 number after decimal point), should be enough
        mxz=num2str(smt(is,3));  % transfer the seismic moment mxz to string, use default precision (4 number after decimal point), should be enough
        myy=num2str(smt(is,4));  % transfer the seismic moment myy to string, use default precision (4 number after decimal point), should be enough
        myz=num2str(smt(is,5));  % transfer the seismic moment myz to string, use default precision (4 number after decimal point), should be enough
        mzz=num2str(smt(is,6));  % transfer the seismic moment mzz to string, use default precision (4 number after decimal point), should be enough
        razi=num2str(srazi(is,ir)); % transfer the source-receiver azimuth to string. For azimuth in degree, using default precision is enough
        ssdp=num2str(soup(is,3),sftp); % transfer this source depth to string
        ssrf=num2str(sroff(is,ir),sftp); % transfer this source-receiver offset to string, no need to add 'blank space' here
        ssyn=sprintf('syn -M%e/%s/%s/%s/%s/%s/%s -A%s -Osdata/S%dR%d.o -Gmodel_%s/%s_%s.grn.0',sm0(is),mxx,mxy,mxz,myy,myz,mzz,razi,is,ir,ssdp,ssrf,srez);
        [~,~]=system(ssyn,'-echo');
    end
end

% use 'readsacp' to trim the synthetic data -- remove the sac data head
for  is=1:nsr
    for ir=1:nre
        stmsacz=sprintf('readsacp file=sdata/S%dR%d.zi et=%d',is,ir,ENt); % note sac need to be installed
        stmsacr=sprintf('readsacp file=sdata/S%dR%d.ri et=%d',is,ir,ENt); % note sac need to be installed
        stmsact=sprintf('readsacp file=sdata/S%dR%d.ti et=%d',is,ir,ENt); % note sac need to be installed
        [~,~]=system(stmsacz,'-echo'); % trim the Z component (Up in fk)
        [~,~]=system(stmsacr,'-echo'); % trim the R component (Radial in fk)
        [~,~]=system(stmsact,'-echo'); % trim the T component (Transverse-clockwise in fk)
    end
end

Nt=ENt; % make the Nt and ENt the same, format the data output
% 1- read the trimmed binary data into matlab
% 2- transfer the ZRT formate to ZNE format (Cylindrical to Cartesian)
% 3- convolve with source time function for each trace
% 4- apply time delay (origin time) for each trace
% 5- sum up over all the sources to obtain the final trace data
wzc=zeros(Nt,nre,nsr); % used to store the Z component data for all sources
wnc=zeros(Nt,nre,nsr); % used to store the N component data for all sources
wec=zeros(Nt,nre,nsr); % used to store the E component data for all sources
fzc=zeros(Nt,nre); % used to store the Z component of the final data
fnc=zeros(Nt,nre); % used to store the N component of the final data
fec=zeros(Nt,nre); % used to store the E component of the final data
for is=1:nsr
    nt0=round(st0(is)/dt)+1; % transfer origin time to number of time sample for each source
    for ir=1:nre
        % read Z component
        filez=sprintf('sdata/S%dR%d.zit',is,ir);
        fidz=fopen(filez,'r');
        zdata=fread(fidz,Nt,'single');
        fclose(fidz);
        % read R component
        filer=sprintf('sdata/S%dR%d.rit',is,ir);
        fidr=fopen(filer,'r');
        rdata=fread(fidr,Nt,'single');
        fclose(fidr);
        % read T component
        filet=sprintf('sdata/S%dR%d.tit',is,ir);
        fidt=fopen(filet,'r');
        tdata=fread(fidt,Nt,'single');
        fclose(fidt);
        % transfer ZRT format to ZNE format
        zc=-zdata; % note in fk the Z component is up, but here we use vertical down
        nc=rdata*cosd(srazi(is,ir))-tdata*sind(srazi(is,ir));
        ec=rdata*sind(srazi(is,ir))+tdata*cosd(srazi(is,ir));
        % convolve with source time function for each trace
        cvzc=conv(zc,stf(is,:)); % Z component
        cvnc=conv(nc,stf(is,:)); % N component
        cvec=conv(ec,stf(is,:)); % E component
        % apply time delay (original time) for each trace
        tdzc=zeros(Nt,1);tdzc(nt0:Nt)=cvzc(1:(Nt-nt0+1)); % Z component
        tdnc=zeros(Nt,1);tdnc(nt0:Nt)=cvnc(1:(Nt-nt0+1)); % N component
        tdec=zeros(Nt,1);tdec(nt0:Nt)=cvec(1:(Nt-nt0+1)); % E component
        % store data into a matrix (dimension 1: time points; dimension 2: all the receivers; dimension 3: all the sources)
        wzc(:,ir,is)=tdzc;
        wnc(:,ir,is)=tdnc;
        wec(:,ir,is)=tdec;
    end
    % sum up over all the sources to obtain the final trace data
    fzc=fzc+wzc(:,:,is);
    fnc=fnc+wnc(:,:,is);
    fec=fec+wec(:,:,is);
end

end
