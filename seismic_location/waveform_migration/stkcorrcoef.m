function scc=stkcorrcoef(data)
% This function is used to calculate the cross-correlation coefficients of
% the input data matrix, and then stack and normalize the correlation coefficients.
%
% INPUT--------------------------------------------------------------------
% data: input data matrix, shape: nt*nre;
%
% OUTPUT-------------------------------------------------------------------
% scc: scalar, the final stacked correlation coefficient;

nre=size(data,2); % get the number of stations

mcc=triu(mycorrcoef(data),1); % correlation coefficient matrix, only keep the up-trianle part, others set to 0

TF = ~isfinite(mcc); % find the index of the NAN/Inf values in the correlation coefficient matrix

mcc(TF)=0; % set the NAN/Inf values to 0

n_nan=sum(TF(:)); % get the number of NAN/Inf values

scc=sum(abs(mcc(:)))/(0.5*nre*(nre-1)-n_nan); % the final normalize stacked CC

end