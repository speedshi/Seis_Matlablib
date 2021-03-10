function cc_sum = runFMF_mdata(components,sname_pick,weights,templates,moveouts,sdata_tag,ofile)
% This funtion is to format the seismic data and run FMF.
% INPUT:
%      components: cell array, specify the seismic data components for
%                  using, e.g.: {'Z','N','E'};
%      sname_pick: cell array, the name of the stations in the templates,
%                  e.g.: {'NRCA'; 'FDMO'; 'CESI'; 'CAMP'; 'TERO'};
%      weights: weight matrix, shape: n_components*n_stations_pick*n_templates;
%      templates: template waveforms, shape:
%                 n_samples_template*n_components*n_stations_pick*n_templates;
%      moveouts: moveout matrix, shape:
%                n_components*n_stations_pick*n_templates;
%      sdata_tag: string, path to continuous seismic data;
%      ofile: string, output filename including path, if [], not save file.
%
% OUTPUT:
%       A '.mat' file which contains FMF results.


n_components = length(components);  % total number of components
n_stations_pick = length(sname_pick);  % total number of stations for template


%% load continuous data of this day
for ic = 1:n_components
    % loop over each component for loading data
    fname_seis = cell(n_stations_pick,1);
    for ist = 1:n_stations_pick
        % set file names of seismic data
        fname_seis{ist} = [sdata_tag sprintf('IV.%s..HH%s.SAC',sname_pick{ist},components{ic})];
    end
    
    seismic_temp = read_seis(fname_seis,components{ic});  % load seismic data
    if components{ic} == 'Z'
        seismic = seismic_temp;
    elseif components{ic} == 'N'
        seismic.ndata = seismic_temp.ndata;
    elseif components{ic} == 'E'
        seismic.edata = seismic_temp.edata;
    end
end
seismic.component = components;


%% construct continuous data for scanning
n_samples_data = size(seismic.zdata,2);
data = zeros(n_samples_data, n_components, n_stations_pick);
for ist = 1:n_stations_pick
    % loop over each station which have picks
    stid = ismember(seismic.name, sname_pick{ist});  % index of this station for seismic dataset
    
    if sum(stid) == 0
        % no data for this station
        data(:,:,ist) = 0;  % assign 0
        weights(:,ist) = 0; % set the weights of this station to 0
    else
        % have data for this station
        for ic = 1:n_components
            if components{ic} == 'Z'
                data(:,ic,ist) = seismic.zdata(stid,:)';
            elseif components{ic} == 'N'
                data(:,ic,ist) = seismic.ndata(stid,:)';
            elseif components{ic} == 'E'
                data(:,ic,ist) = seismic.edata(stid,:)';
            end
        end
        
    end
    
end


%% run template matching-----------------------------------------------
step = 1;
weights = weights / sum(weights(:));  % sum(weights) should be 1
[cc_sum] = fast_matched_filter(templates, moveouts, weights, data, step);


%% save results of this day--------------------------------------------
if ofile
    para.t0 = seismic.t0;
    para.dt = seismic.dt;
    save(ofile,'cc_sum','weights','para','-v7.3');
end


end