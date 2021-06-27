function gene_FMFtemplate(Nobs_file, output_dir, seisd_dir, dpara, ops)
%%% This function is used to generate the templates for FMF (fast matched
%%% filter).
% INPUTS:
% Nobs_file: str, filename of the NonLinLoc observation file containing
% phase picks;
% output_dir: str, directory for outputs;
% seisd_dir: str, directory to seismic data;
% dpara: structure, parameters for generating templates, examples:
%       dpara.components = {'Z','N','E'};  % the components of seismic data
%       templates are extract around P-wave on Z component and around S-wave on E/N components
%       dpara.twindl = 2;  % window length in second for templates
%       dpara.windb_p = 0.2;  % window length in second before the P-wave arrival-time
%       dpara.windb_s = 0.2;  % window length in second before the S-wave arrival-time
% ops: structure, other parameters required such as for plotting waveforms.
%       examples:
%       ops.wind_P1 = 1;  % left window length in second, relative P-wave arrival-time
%       ops.wind_P2 = 20;  % right window length in second, relative P-wave arrival-time
%       ops.wind_S1 = 6;  % left window length in second, relative S-wave arrival-time
%       ops.wind_S2 = 15;  % right window length in second, relative S-wave arrival-time
%
% OUTPUT:
% 'template_dataset.mat': dateset containing templates and related
% parameters;
% 'station-name_3C_template_waveforms.png': waveform plots for the
% templates at different stations.


%% get P- and S-phase arrival-times from NLLoc observation data file
pick = getpicks_fromNLLobs(Nobs_file);


%% the picking info--------------------------------------------------------
sname_pick = fieldnames(pick);  % the name of the picked stations
n_stations_pick = length(sname_pick);  % total number of stations for template
n_components = length(dpara.components);  % total number of components


%% prepare templates-------------------------------------------------------
% load three component data
temp1 = struct2cell(pick.(sname_pick{1}));  % get a picking time
srftime = temp1{1};  % reference time for extracting seismic data
sdata_tag = [seisd_dir sprintf('/%d/continuous%03d',year(srftime),day(srftime,'dayofyear')) '/processed_5-25Hz/'];  % path to the waveforms to extract templates
for ic = 1:n_components
    % loop over each component for loading data
    
    fname_seis = cell(n_stations_pick,1);
    for ist = 1:n_stations_pick
        % file names of seismic data
        fname_seis{ist} = [sdata_tag sprintf('IV.%s..HH%s.SAC',sname_pick{ist},dpara.components{ic})];
    end
    seismic_temp = read_seis(fname_seis,dpara.components{ic});  % load seismic data
    if dpara.components{ic} == 'Z'
        seismic = seismic_temp;
    elseif dpara.components{ic} == 'N'
        seismic.ndata = seismic_temp.ndata;
    elseif dpara.components{ic} == 'E'
        seismic.edata = seismic_temp.edata;
    end
end
seismic.component = dpara.components;
dpara.dt = seismic.dt;  % time sampling interval in second
ttmax = seismic.t0 + seconds((seismic.NPTS-1)*seismic.dt);  % ending time of seismic trace

% output figures for showing the waveforms of templates
for ist = 1:n_stations_pick
    % loop over each station which have picks
    try
        tp = getfield(pick,sname_pick{ist},'P');  % get the P-wave arrival-time for a station
    catch
        warning('No P picks for station: %s.\n',sname_pick{ist});
        tp = [];
    end
    try
        ts = getfield(pick,sname_pick{ist},'S');  % get the S-wave arrival-time for a station
    catch
        warning('No S picks for station: %s.\n',sname_pick{ist});
        ts = [];
    end
    if ~isempty(tp)
        timerg = [tp-seconds(ops.wind_P1), tp+seconds(ops.wind_P2)];  % time range for extract the data
    else
        timerg = [ts-seconds(ops.wind_S1), ts+seconds(ops.wind_S2)];  % time range for extract the data
    end
    
    % make sure time range donot exceed data time range
    if timerg(1) < seismic.t0
        warning('Starting time window not appropriate, reset to: %s.\n', seismic.t0);
        timerg(1) = seismic.t0;
    end
    if timerg(2) > ttmax
        warning('Ending time window not appropriate, reset to: %s.\n', ttmax);
        timerg(2) = ttmax;
    end
    
    stid = ismember(seismic.name, sname_pick{ist});  % index of this station for seismic dataset
    fig = figure('visible', 'off');
    for ic = 1:n_components
        if dpara.components{ic} == 'Z'
            sdata = seismic.zdata(stid,:)';
            [sdata, taxis] = seisext(sdata,seismic.dt,seismic.t0,timerg);
            subplot(n_components,1,ic);
            plot(taxis,sdata,'k','linewidth',0.5); hold on;
            if ~isempty(tp)
                xline(tp,'b','linewidth',0.5); hold on;
            end
            if ~isempty(ts)
                xline(ts,'r','linewidth',0.5); hold on;
            end
            ylabel('Z','FontWeight','bold');
            set(gca, 'XLimSpec', 'Tight');
        elseif dpara.components{ic} == 'N'
            sdata = seismic.ndata(stid,:)';
            [sdata, taxis] = seisext(sdata,seismic.dt,seismic.t0,timerg);
            subplot(n_components,1,ic);
            plot(taxis,sdata,'k','linewidth',0.5); hold on;
            if ~isempty(tp)
                xline(tp,'b','linewidth',0.5); hold on;
            end
            if ~isempty(ts)
                xline(ts,'r','linewidth',0.5); hold on;
            end
            ylabel('N','FontWeight','bold');
            set(gca, 'XLimSpec', 'Tight');
        elseif dpara.components{ic} == 'E'
            sdata = seismic.edata(stid,:)';
            [sdata, taxis] = seisext(sdata,seismic.dt,seismic.t0,timerg);
            subplot(n_components,1,ic);
            plot(taxis,sdata,'k','linewidth',0.5); hold on;
            if ~isempty(tp)
                xline(tp,'b','linewidth',0.5); hold on;
            end
            if ~isempty(ts)
                xline(ts,'r','linewidth',0.5); hold on;
            end
            ylabel('E','FontWeight','bold');
            set(gca, 'XLimSpec', 'Tight');
        end
    end
    sgtitle(['Station: ' sname_pick{ist}],'FontWeight','bold');
    figname = sprintf('%s/%s_3C_template_waveforms',output_dir,sname_pick{ist});
    print(figname,'-dpng','-r300');
    close(fig);
end

% construct templates
% templates are extract around P-wave on Z component and around S-wave on E/N components
n_samples_template = round(dpara.twindl / seismic.dt) + 1;
n_templates = 1;
templates = zeros(n_samples_template,n_components,n_stations_pick,n_templates);  % templates
moveouts = NaN(n_components, n_stations_pick, n_templates);  % moveouts
weights = ones(n_components, n_stations_pick, n_templates);  % weights
for ist = 1:n_stations_pick
    % loop over each station which have picks
    stid = ismember(seismic.name, sname_pick{ist});  % index of this station for seismic dataset
    
    try
        tp = getfield(pick,sname_pick{ist},'P');  % get the P-wave arrival-time for a station
    catch
        tp = [];
    end
    try
        ts = getfield(pick,sname_pick{ist},'S');  % get the S-wave arrival-time for a station
    catch
        ts = [];
    end
    
    if ~isempty(tp)
        tp_id = round(seconds(tp-seconds(dpara.windb_p)-seismic.t0)/seismic.dt)+1;  % index for P-wave window starting time
        twrgid_p = [tp_id, tp_id+n_samples_template-1];  % sample index range for extracting data around P-wave
    else 
        tp_id = [];
        twrgid_p = [];
    end
    if ~isempty(ts)
        ts_id = round(seconds(ts-seconds(dpara.windb_s)-seismic.t0)/seismic.dt)+1;  % index for S-wave window starting time
        twrgid_s = [ts_id, ts_id+n_samples_template-1];  % sample index range for extracting data around S-wave
    else
        ts_id = [];
        twrgid_s = [];
    end
    
    % make sure date extracting time range donot exceed limit
    if ~isempty(twrgid_p) && (twrgid_p(1) < 1)
        warning('Starting time for extracting P-phase preceeds trace t0! Neglect P-phase for stations: %s!\n', sname_pick{ist});
        tp_id = [];
        twrgid_p = [];
    end
    if ~isempty(twrgid_p) && (twrgid_p(2) > seismic.NPTS)
        warning('Ending time for extracting P-phase exceeds trace ending time! Neglect P-phase for stations: %s!\n', sname_pick{ist});
        tp_id = [];
        twrgid_p = [];
    end
    if ~isempty(twrgid_s) && (twrgid_s(1) < 1)
        warning('Starting time for extracting S-phase preceeds trace t0! Neglect S-phase for stations: %s!\n', sname_pick{ist});
        ts_id = [];
        twrgid_s = [];
    end
    if ~isempty(twrgid_s) && (twrgid_s(2) > seismic.NPTS)
        warning('Ending time for extracting S-phase exceeds trace ending time! Neglect S-phase for stations: %s!\n', sname_pick{ist});
        ts_id = [];
        twrgid_s = [];
    end
    
    for ic = 1:n_components
        if dpara.components{ic} == 'Z'
            if ~isempty(tp_id)
                templates(:,ic,ist,1) = seismic.zdata(stid,twrgid_p(1):twrgid_p(2))';
                moveouts(ic,ist,1) = tp_id;  % moveouts for P-waves
            else
                % no picking for P-phase
                templates(:,ic,ist,1) = zeros(n_samples_template,1);
                moveouts(ic,ist,1) = NaN;  % moveouts for P-waves, cannot assign 0 here as we need to find the minimal later
                weights(ic,ist,1) = 0;  % weights for P-waves
            end
        elseif dpara.components{ic} == 'N'
            if ~isempty(ts_id)
                templates(:,ic,ist,1) = seismic.ndata(stid,twrgid_s(1):twrgid_s(2))';
                moveouts(ic,ist,1) = ts_id;  % moveouts for S-waves
            else
                % no picking for S-phase
                templates(:,ic,ist,1) = zeros(n_samples_template,1);
                moveouts(ic,ist,1) = NaN;
                weights(ic,ist,1) = 0;
            end
        elseif dpara.components{ic} == 'E'
            if ~isempty(ts_id)
                templates(:,ic,ist,1) = seismic.edata(stid,twrgid_s(1):twrgid_s(2))';
                moveouts(ic,ist,1) = ts_id;  % moveouts for S-waves
            else
                % no picking for S-phase
                templates(:,ic,ist,1) = zeros(n_samples_template,1);
                moveouts(ic,ist,1) = NaN;
                weights(ic,ist,1) = 0;
            end
        end
    end
end
moveouts_original = moveouts;  % save the original moveout data
moveouts = moveouts - min(moveouts(:));  % let the minimal moveout to be 0
moveouts(isnan(moveouts)) = 0;  % replace NaN with 0


%% save templates and the related info-------------------------------------
save([output_dir '/template_dataset.mat'], 'templates', 'moveouts', 'moveouts_original', 'weights', 'pick', 'dpara','-v7.3');


end
