function wavefmsk(fzc,travelp,travels,st0r,t0dt,nsxr,nsyr,nszr,dt)
% This function is used to performing waveform stacking.
% All the input and output data are stored and saved at folder '../data'.


% calculate and set some parameters
nre=size(travelp,2); % number of geophones in travel-time calculation
Nt=size(fzc,1); % number of time samples
st0=st0r(1):t0dt:st0r(2); % original time range
nst0=max(size(st0)); % number of searching points

for ik=6:7
    % transfer the full waveform data using some stack kernel functions
    if ik==0
        % use original data;
        skfz=stack_kernelf(fzc,0);
    elseif ik==1
        % use absolute data;
        skfz=stack_kernelf(fzc,1);    
    elseif ik==3
        % use envelope data;
        skfz=stack_kernelf(fzc,2);
    elseif ik==4
        % use non-negtive data;
        skfz=stack_kernelf(fzc,3);
    elseif ik==5
        % use STA/LTA
        skfz=my_stalta(fzc,10,100,0,0); % square data - classical STA/LTA
    elseif ik==6
        % use Kurtosis
        skfz=my_kurtosis(fzc,100); % calculate Kurtosis within a perticular time window
    elseif ik==7
        % use positive Kurtosis derivative
        ktoz=my_kurtosis(fzc,100); % calculate Kurtosis within a perticular time window
        skfz=intder(ktoz,dt,0); % calculate the derivatives of the Kurtosis
        skfz(skfz<0)=0; % set the negative derivatives to be 0
    else
        error('Wrong input for stack kernel function!');
    end
    
    % stack the data and obtain the imaging results
    wfmstk=zeros(nst0,nsxr,nsyr,nszr); % array used to save the stacked data, demision 1: time; demision 2: X; demision 3: Y; demision 4: Z
    for iz=1:nszr
        for iy=1:nsyr
            for ix=1:nsxr
                for it=1:nst0
                    id=(iz-1)*nsxr*nsyr+(iy-1)*nsxr+ix; % index for finding the source in 'tvtsoup' and travel-time in 'travelp'
                    tvpn=round((travelp(id,:)+st0(it))/dt)+1; % time point of direct P-wave for this source position
                    tvsn=round((travels(id,:)+st0(it))/dt)+1; % time point of direct S-wave for this source position
                    for ir=1:nre
                        if  (tvpn(ir)<=Nt && tvpn(ir)>=1)
                            fmp=skfz(tvpn(ir),ir);
                        else
                            fmp=0;
                        end
                        if (tvsn(ir)<=Nt && tvsn(ir)>=1)
                            fms=skfz(tvsn(ir),ir);
                        else
                            fms=0;
                        end
                        wfmstk(it,ix,iy,iz)=wfmstk(it,ix,iy,iz)+fmp+fms;
                    end
                    
                end
            end
        end
    end
    
    % save the results
    dname=['imagz' num2str(ik) '_ns0d5h']; % data name
    eval([dname '=wfmstk;']); % rename the variable
    eval(['save ../data/homo/' dname '.mat ' dname ' -v7.3;']) % save the data
end

end
