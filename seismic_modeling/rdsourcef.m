function [dt, nt, ns, source]=rdsourcef(fname)
% Read the source file.
% Note the input moment tensor will be normalized, which means the square
% root of the sum of the squares of the elements of each of the moment
% tensor is equal to one.
% INPUT:-------------------------------------------------
% fname: file name of the source file.
% OUTPUT:----------------------------------------------
% dt: time sampling interval of the modeling and source time function;
% nt: number of time samples of the modeling and source time function;
% ns: number of sources;
% source: structure array containing source information;
% source.pos: X-Y-Z positions of source, 3*1;
% source.stf: source time function of source, nt*1;
% source.sftp: source time function type, scalar;
% source.freq: dominant frequency of the source time function, scalar;
% source.t0: origin time or delay time of the source, scalar;
% source.smt: moment tensor of source, 3*3;
% source.m0: scalar seismic moment, unit: N*m.


if nargin<1
    fname="source.dat";
end

fid1=fopen(fname,'r');

cdt=textscan(fid1,"%f",1);
dt=cdt{1}; % time sampling interval of the modeling and source time function

cnt=textscan(fid1,"%d",1);
nt=double(cnt{1}); % number of time samples of the modeling and source time function

snum=textscan(fid1,"%d",1);
ns=double(snum{1}); % number of sources

% preallocate a structure array
source=struct('pos',zeros(3,1),'stf',zeros(nt,1),'sftp',zeros(1,1),'freq',zeros(1,1),'t0',zeros(1,1),'smt',zeros(3,3),'m0',zeros(1,1));

for ii=1:ns
    sloc=textscan(fid1,"%f %f %f",1); % source location
    source(ii).pos(1,1)=sloc{1,1}; % X
    source(ii).pos(2,1)=sloc{1,2}; % Y
    source(ii).pos(3,1)=sloc{1,3}; % Z
    stftp=textscan(fid1,"%d",1); % type of the source time function
    source(ii).sftp=double(stftp{1}); % source time function type
    if stftp{1}==0
        % read wavelet from file
        wvlti=textscan(fid1,"%s",1); % the name of the wavelet file
        ogtime=textscan(fid1,'%f',1); % origin time of the source
        source(ii).t0=ogtime{1}; % origin time or delayed time
        source(ii).freq=nan; % frequency of source time function, not known
        stemp=load(wvlti{1}{1}); % read in source time function from file
        nn=size(stemp(:),1);
        n0=round(ogtime{1}/dt); % number of samples used to delay the source time function
        if n0+nn<nt
            source(ii).stf(1:n0,1)=0;
            source(ii).stf(n0+1:n0+nn,1)=stemp;
            source(ii).stf(nn+1:nt,1)=0;
        elseif n0+nn==nt
            source(ii).stf(1:n0,1)=0;
            source(ii).stf(n0+1:nt,1)=stemp;
        else
            source(ii).stf(1:n0,1)=0;
            source(ii).stf(n0+1:nt,1)=stemp(1:nt-n0);
            warning(['The source time function of source %d is too long! ' ...
                'The total number of time samples in this modeling is %d. ' ...
                'However the source time funciton of source %d has %d time samples + %d origin time sample shifts. ' ...
                'Transection of the source time function has to be made!'],ii,nt,ii,nn,n0);
        end
    else
        % theoretical wavelet
        wvlti=textscan(fid1,"%f",1); % dominant frequency
        ogtime=textscan(fid1,'%f',1); % origin time of the source
        source(ii).t0=ogtime{1}; % origin time or delayed time
        source(ii).freq=wvlti{1}; % frequency of source time function
        if stftp{1}==1
            % use Ricker wavelet
            source(ii).stf=rickerw(wvlti{1},dt,nt,ogtime{1});
        elseif stftp{1}==-1
            % use impulse as source time function, calculate the Green's function
            source(ii).freq=nan; % reset the frequency
            source(ii).stf=zeros(size(nt,1));
            n0=round(ogtime{1}/dt); % number of samples used to delay the source time function
            if n0>0
                source(ii).stf(n0,1)=1; % set impulse time
            else
                source(ii).stf(1,1)=1; % set impulse time
            end
        end
        
    end
    sclmt=textscan(fid1,"%f",1); % scalar moment
    source(ii).m0=sclmt{1}; % scalar seismic moment, in N*m
    itpmt=textscan(fid1,"%d",1);% the input type of moment tensor
    if itpmt{1}==0
        % the use of fault geometry parameters: strike, rake and dip
        mto=textscan(fid1,"%f %f %f",1);
        source(ii).smt=fgeom2mt(mto{1},mto{2},mto{3});
    else
        % directly input moment tensor
        mto=textscan(fid1,"%f %f %f %f %f %f",1);
        source(ii).smt(1,1)=mto{1}; % mxx
        source(ii).smt(2,2)=mto{2}; % myy
        source(ii).smt(3,3)=mto{3}; % mzz
        source(ii).smt(1,2)=mto{4}; % mxy
        source(ii).smt(1,3)=mto{5}; % mxz
        source(ii).smt(2,3)=mto{6}; % myz
        source(ii).smt(2,1)=source(ii).smt(1,2); % myx
        source(ii).smt(3,1)=source(ii).smt(1,3); % mzx
        source(ii).smt(3,2)=source(ii).smt(2,3); % myz
    end
    % normalize the input seismic moment tensor
    source(ii).smt=mtnorm(source(ii).smt);
end

fclose(fid1);

end