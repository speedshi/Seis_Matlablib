function [ux,uy,uz,dt]=gsynwihl(sname,rname,mname)
% This function is to generate synthetic waveforms in a homogeneous or
% layered model. The output waveforms are the particle displacements in meter.
% INPUT:-------------------------------------------------
% sname: the file name of source file;
% rname: the file name of receiver file;
% mname: the file name of velocity model file.
% OUTPUT:----------------------------------------------
% ux: displacement in the X direction, nt*nr, in meter;
% uy: displacement in the Y direction, nt*nr, in meter;
% uz: displacement in the Z direction, nt*nr, in meter;
% dt: time sampling interval of the output traces, in second.


if nargin < 1
    sname="source.dat";
    rname="receiver.dat";
    mname="model.dat";
end

% read the source file
[dt,nt,ns,source]=rdsourcef(sname);

% read the receiver file
[nr,recp]=rdreceiverf(rname);

% read the velocity model file
[nl,rsd0,model]=rdmodelf(mname);

ux=zeros(nt,nr); % X component of wavefields
uy=zeros(nt,nr); % Y component of wavefields
uz=zeros(nt,nr); % Z component of wavefields

if nl==1
    % homogeneous model
    for ids=1:ns
        % generate wavefields of the 'ids'-th source
        
        % integrate the scalar seismic moment into the moment tensor
        smt=source(ids).smt*source(ids).m0;
    
        if source(ids).sftp==1
            % use Ricker wavelet, can generate waveforms using the
            % analytical expressions of the wavelet
            [tux,tuy,tuz]=gsynwhomo_rickerw(model.vp,model.vs,model.den,recp,source(ids).pos,smt,dt,nt,source(ids).freq);
        elseif source(ids).sftp==-1
            % use impuse as source time function, calculate the Green's function
            [tux,tuy,tuz]=homogreenf(model.vp,model.vs,model.den,recp,source(ids).pos,smt,dt,nt,source(ids).t0);
        else
            % use other wavelet or input source time function
            [tux,tuy,tuz]=gsynwhomo(model.vp,model.vs,model.den,recp,source(ids).pos,smt,dt,nt,source(ids).stf);
        end
        
        % add the wavefields of the source ids to the whole wavefield
        ux=ux+tux;
        uy=uy+tuy;
        uz=uz+tuz;
    end
else
    % layered model, use FK to generate waveforms
    
    % check and make sure the depths of sources and receivers are within
    % the model definition
    sdept=zeros(ns,1); % source depths
    for ids=1:ns
        sdept(ids)=source(ids).pos(3);
    end
    minsdep=min(sdept); % the minimal source depth, should be larger than the starting depth of the velocity model
    minrdep=min(recp(:,3)); % the minimal receiver depth, should be larger than the starting depth of the velocity model
    if minsdep<rsd0
        warning('Source depth out of range. Minimal source depth (%f m) is shallower than the input starting depth (%f m) of the velocity model.',minsdep,rsd0);
    end
    if minrdep<rsd0
        warning('Receiver depth out of range. Minimal receiver depth (%f m) is shallower than the input starting depth (%f m) of the velocity model.',minrdep,rsd0);
    end
    mindep=min([minsdep minrdep]); % the minimal source and receiver depth
    if mindep<rsd0
        thick1=model.thickness(1);
        stdpo=rsd0;
        model.thickness(1)=model.thickness(1)+rsd0-mindep; % adjust the thickness of the first layer to include the minimal depth
        rsd0=mindep; % reset the starting depth (free surface) to the minimal depth
        warning('In order to include all sources and receivers in the velocity model, the thickness of the first layer is reset from %f to %f, the starting depth is reset from %f to %f.',thick1,model.thickness(1),stdpo,rsd0);
    end
    
    % generate model input file for FK
    wrtmdf(model.thickness,model.vp,model.vs,model.den,model.qp,model.qs); % note the units of the input and output parameters
    soupk=zeros(ns,3); % source locations
    recpk=recp/1000; % receiver locations, unite: m->km.
    smt=zeros(ns,6); % source moment tensors, note the order of moment tensor components
    sstf=zeros(ns,nt); % source time functions
    sm0=zeros(ns,1); % scalar seismic moment
    utran=1e7; % used to tranfer the unit of moment tensor from N*m to dyne*cm     
    for ids=1:ns
        soupk(ids,:)=source(ids).pos/1000; % source location, X-Y-Z, note unit transformation, m->km.
        % Note for the input parameter file 'source.dat', the unit of moment tensor is N*m; but for fk, it uses dyne*cm. So we need to transfer the unit.
        sm0(ids,1)=source(ids).m0*utran; % scalar seismic moment, note the unit transfer
        smt(ids,1)=source(ids).smt(1,1); % mxx
        smt(ids,2)=source(ids).smt(1,2); % mxy
        smt(ids,3)=source(ids).smt(1,3); % mxz
        smt(ids,4)=source(ids).smt(2,2); % myy
        smt(ids,5)=source(ids).smt(2,3); % myz
        smt(ids,6)=source(ids).smt(3,3); % mzz
        sstf(ids,:)=source(ids).stf; % source time functions
    end
    % Note the fk assume the depth of the free surface is 0. So here we correct for the depths of sources and receivers according to the input starting depth.
    soupk(:,3)=soupk(:,3)-rsd0/1000; % correct for the source depths. Note unit transfer.
    recpk(:,3)=recpk(:,3)-rsd0/1000; % correct for the receiver depths. Note unit transfer.
    
    [~,~,~,uz,ux,uy]=gsynwavefk(soupk,recpk,smt,sm0,sstf,dt,nt); % call fk to generate wavefields
    % For fk, the output is displacement in cm. Here we transform the unit
    % of displacement from cm to m.
    tudis=0.01; % used to transfer the unit from cm to m
    ux=ux*tudis;
    uy=uy*tudis;
    uz=uz*tudis;
end


end