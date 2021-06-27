function pick = getpicks_fromNLLobs(Nobs_file)
% This function is used to obtain picks (P- and S-phases arrivaltimes) from
% NonLinLoc observation file.
% INPUT:
%       Nobs_file: filename of the input NonLinLoc obervation file. Check
%       NonLinLoc software for the format of this file.
% OUTPUT:
%       pick: structure, contains the picked arrivaltimes of P- and/or
%       S-phases at different stations;
%       example: 
%               pick.stationname.P: P-phase arrivaltime at station: 'stationame';  
%               pick.stationname.S: S-phase arrivaltime at station: 'stationame'.


fid1 = fopen(Nobs_file,'r');
dd = textscan(fid1,'%s %*s %*s %*s %s %*s %s %s %s %*s %*f %*f %*f %*f');
fclose(fid1);

nphs = size(dd{1},1);  % the total number of phases picked
for ii = 1:nphs
    pick.(dd{1}{ii}).(dd{2}{ii}) = datetime([dd{3}{ii} dd{4}{ii} dd{5}{ii}],'InputFormat','yyyyMMddHHmmss.SSSS');
end


end
