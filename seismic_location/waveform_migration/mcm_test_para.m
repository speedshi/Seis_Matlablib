function [s_pro,n_var]=mcm_test_para(trace,mcm,search,earthquake,fshow)
% This function is used to perform tests and decide MCM parameters.
%
% INPUT--------------------------------------------------------------------
% trace: matlab structure, contains seismic data information;
% trace.data: seismic data, 2D array, nre*nt;
% trace.dt: time sampling interval, in second;
% trace.travelp: P-wave traveltime table, 2D array, nsr*nre;
% trace.travels: S-wave traveltime table, 2D array, nsr*nre;
% mcm: matlab structure, contains MCM parameters;
% mcm.phasetp: specify seismic phase used for migration, scalar;
% mcm.tpwind: P-phase time window length in second, scalar;
% mcm.tswind: S-phase time window length in second, scalar;
% search: matlab structure, describe the imaging area,
% search.soup: source imaging positions, 2D array, ns*3, in meter;
% earthquake: matlab structure, contains the location and origin time of the earthquake;
% earthquake.north: scalar, earthquake location in north direction, in meter;
% earthquake.east: scalar, earthquake location in east direction, in meter;
% earthquake.depth: scalar, earthquake location in depth direction, in meter;
% earthquake.t0: scalar, relative earthquake origin time, in second,
% relative to the origin time of the seismic data;
% fshow: logical, whether to plot the figures, default is to plot;
%
% OUTPUT-------------------------------------------------------------------
% s_pro: source prominence;
% n_var: noise variance;

% set default value
if nargin<5
    fshow=true;
end

% the location of the earthquake, in meter
soup_cata=[earthquake.north earthquake.east earthquake.depth];

% obtain the origin time of the tested event relative to the starting time of traces (t0), in second
stime=earthquake.t0;

% determine catalog event index in the soup, by using minimal distance--approximate
[~,idseca]=min(sum((search.soup-soup_cata).^2,2));

switch mcm.phasetp
    case 0
        % only for P-waves
        fprintf('Extract and align the P-waves for different traces.\n');
        twinl=90; % left window size, in second
        twinr=300; % right window seize, in second
        
        % build arrival times for the tested event, in second
        atimes=stime+trace.travelp(idseca,:);
        
        % extrace the data
        edata=wave_extract(trace.data',trace.dt,atimes,twinl,twinr);
        
        if fshow
            % display the extracted data, note 'edata' has the shape of 'nt*nre'
            seisrsdisp(edata,trace.dt);
            hold on; ax1=gca; ylim1=ax1.YLim;
            plot([twinl,twinl],ylim1,'r'); hold on;
        end
        
    case 1
        % only for S-waves
        fprintf('Extract and align the S-waves for different traces.\n');
        twinl=90; % left window size, in second
        twinr=300; % right window seize, in second
        
        % build arrival times for the tested event, in second
        atimes=stime+trace.travels(idseca,:);
        
        % extrace the data
        edata=wave_extract(trace.data',trace.dt,atimes,twinl,twinr);
        
        if fshow
            % display the extracted data, note 'edata' has the shape of 'nt*nre'
            seisrsdisp(edata,trace.dt);
            hold on; ax1=gca; ylim1=ax1.YLim;
            plot([twinl,twinl],ylim1,'r'); hold on;
        end
        
    case 2
        % for both P- and S-waves
        % the P- and S-waves are extract separately and then concatenate
        % together
        fprintf('Extract and align both the P and S-waves for different traces.\n');
        
        % for P-wave segment
        twinl1=90; % left window size, in second
        twinr1=3; % right window seize, in second
        
        % build arrival times for the tested event, in second
        atimes=stime+trace.travelp(idseca,:);
        
        % extrace the data
        edata1=wave_extract(trace.data',trace.dt,atimes,twinl1,twinr1);
        
        
        % for S-wave segment
        twinl2=0.5; % left window size, in second
        twinr2=300; % right window seize, in second
        
        % build arrival times for the tested event, in second
        atimes=stime+trace.travels(idseca,:);
        
        % extrace the data
        edata2=wave_extract(trace.data',trace.dt,atimes,twinl2,twinr2);
        
        % combine P and S segments together
        edata=[edata1; edata2];
        
        if fshow
            % display the extracted data, note 'edata' has the shape of 'nt*nre'
            seisrsdisp(edata,trace.dt);
            hold on; ax1=gca; ylim1=ax1.YLim;
            plot([twinl1,twinl1],ylim1,'r'); hold on;
            stpt=twinl1+twinr1+twinl2;
            plot([stpt,stpt],ylim1,'r'); hold on;
        end
        
end

% show the seismogram and spectrogram
i_station=1; % specify to show which station
ispectrogram_1(edata(:,i_station),trace.dt,trace.name{i_station},2);

[nt,~]=size(edata); % obtain the number of time samples and stations

switch mcm.phasetp
    case 0
        % for P-phase
        twind=mcm.tpwind; % migration time window
        n_win=round(twind/trace.dt); % number of time points in the migration window
        
        nn=nt-n_win+1; % number of migration time points
        migv=zeros(nn,1); % initialize the migration value array
        
        for it=1:nn
            migv(it)=stkcorrcoef(edata(it:it+n_win-1,:));
        end
        
        if fshow
            % show the migration values
            figure;
            plot((0:nn-1)*trace.dt,migv,'k','linewidth',1.5); hold on;
            ax=gca; ylim=ax.YLim;
            plot([twinl twinl],ylim,'r'); hold on;% origin time (arrival time of P)
            t_cor=twinl-twind+0.5; % calibrate according to time window length and period
            plot([t_cor t_cor],ylim,'r--'); hold on;% calibrated origin time (arrival time of P)
            xlabel('Time (s)'); ylabel('Coherency');
            title('Migration trace using P-phase'); axis tight;
        end
        
    case 1
        % for S-phase
        twind=mcm.tswind; % migration time window
        n_win=round(twind/trace.dt); % number of time points in the migration window
        
        nn=nt-n_win+1; % number of migration time points
        migv=zeros(nn,1); % initialize the migration value array
        
        for it=1:nn
            migv(it)=stkcorrcoef(edata(it:it+n_win-1,:));
        end
        
        if fshow
            % show the migration values
            figure;
            plot((0:nn-1)*trace.dt,migv,'k','linewidth',1.5); hold on;
            ax=gca; ylim=ax.YLim;
            plot([twinl twinl],ylim,'r'); hold on;% origin time (arrival time of S)
            t_cor=twinl-twind+0.5; % calibrate according to time window length and period
            plot([t_cor t_cor],ylim,'r--'); hold on;% calibrated origin time (arrival time of S)
            xlabel('Time (s)'); ylabel('Coherency');
            title('Migration trace using S-phase'); axis tight;
        end
        
    case 2
        % for P- and S-phase
        twind=mcm.tpwind; % migration time window
        n_winp=round(twind/trace.dt); % number of time points in the migration window for P
        
        twind2=mcm.tswind; % migration time window
        n_wins=round(twind2/trace.dt); % number of time points in the migration window for S
        
        n_intv=round((twinr1+twinl2)/trace.dt); % number of time points between the P- and S-phase
        
        nn=nt-n_intv-n_wins; % number of migration time points
        migv=zeros(nn,1); % initialize the migration value array
        
        for it=1:nn
            cc_p=stkcorrcoef(edata(it:it+n_winp-1,:)); % stacked CC of P
            cc_s=stkcorrcoef(edata((it+n_intv+1):(it+n_intv+1+n_wins-1),:)); % stacked CC of S
            migv(it)=0.5*(cc_p+cc_s); % final migration value
        end
        
        if fshow
            % show the migration values
            figure;
            plot((0:nn-1)*trace.dt,migv,'k','linewidth',1.5); hold on;
            ax=gca; ylim=ax.YLim;
            plot([twinl1 twinl1],ylim,'r'); hold on;% origin time (also P arrival time)
            t_cor=twinl1-0.5*(mcm.tpwind+mcm.tswind)+0.5; % calibrate according to time window length and period
            plot([t_cor t_cor],ylim,'r--'); hold on;% calibrated origin time
            plot([twinl1+twinr1+twinl2 twinl1+twinr1+twinl2],ylim,'b'); hold on;% S arrival time
            xlabel('Time (s)'); ylabel('Coherency');
            title('Migration trace using P- and S-phases'); axis tight;
        end
        
end

[mm,mm_id]=max(migv); % obtain the maximum migration value and time

if fshow
    % show maximum coherent time window on record sections
    wst=(mm_id-1)*trace.dt;
    xx=[wst wst+twind wst+twind wst];
    yy=[ylim1(1) ylim1(1) ylim1(2) ylim1(2)];
    fill(ax1,xx,yy,'y','linestyle','none','FaceAlpha',0.35); hold on;
    
    if mcm.phasetp==2
        wst=(mm_id-1)*trace.dt+twinr1+twinl2;
        xx=[wst wst+twind2 wst+twind2 wst];
        fill(ax1,xx,yy,'g','linestyle','none','FaceAlpha',0.35); hold on;
    end
end

% calculate source prominence
s_pro=mm/median(migv);

% calculate noise variance
if mcm.phasetp==2
    sdd=round((twinl1+20)/trace.dt);
else
    sdd=round((twinl+20)/trace.dt);
end
n_var=var(migv(sdd:end));

if fshow
    fprintf('Source prominence: %f.\n',s_pro);
    fprintf('Noise variance: %f.\n',n_var);
end

end