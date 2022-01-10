% MSAC toolkit function.
%
%  [sac_trace_structure] = msac_read(filename) ;
%
%  Read a binary-formatted sac file, and return a structure containing the
%  header elements, and data
%

%-------------------------------------------------------------------------------
%
%  This software is distributed under the term of the BSD free software license.
%
%  Copyright:
%     (c) 2003-2010, James Wookey
%
%  All rights reserved.
%
%   * Redistribution and use in source and binary forms, with or without
%     modification, are permitted provided that the following conditions are
%     met:
%        
%   * Redistributions of source code must retain the above copyright notice,
%     this list of conditions and the following disclaimer.
%        
%   * Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the
%     documentation and/or other materials provided with the distribution.
%     
%   * Neither the name of the copyright holder nor the names of its
%     contributors may be used to endorse or promote products derived from
%     this software without specific prior written permission.
%
%
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
%   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
%   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
%   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
%   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
%   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
%   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%-------------------------------------------------------------------------------
 
function [sactrace] = msac_read(fname) ;

[fid,message] = fopen(fname,'r') ;

if (~isempty(message))
   fprintf('File: %s\n',fname) ;
   fprintf('%s\n',message) ;
   return ;
end

% get the float header
hd = fread(fid,70,'float32') ;

% get the integer header
ihd = fread(fid,40,'int') ;

[dum,dum2,endian_str] = computer ;

% Test for wrong endian-format, using version number in header, should be 6
% but it will be something bizarre in the wrong endian format
if (ihd(007)<1 | ihd(007)>10)
%  Reread with reversed endian format
   fclose(fid) ;

% get computer type (big or little endian)
% reopen file with reversed endian set
   if ( endian_str == 'L' )
      fid = fopen(fname,'r','b') ;
      tmpEndian = 'b' ;
   else
      fid = fopen(fname,'r','l') ;
      tmpEndian = 'l' ;
   end  
   
   % get the float header
   hd = fread(fid,70,'float32') ;

   % get the integer header
   ihd = fread(fid,40,'int') ;
   
else
   tmpEndian = lower(endian_str) ;

end

% check that the header version number is valid
if (ihd(007) ~= 6) 
   error('SAC file %s header version is out of date; ie., not equal to 6',fname) ;
end

% convert the real header

sactrace.delta     = hd(001) ;
sactrace.depmin    = hd(002) ;
sactrace.depmax    = hd(003) ;
sactrace.scale     = hd(004) ;
sactrace.odelta    = hd(005) ;
sactrace.b         = hd(006) ;
sactrace.e         = hd(007) ;
sactrace.o         = hd(008) ;
sactrace.a         = hd(009) ;
sactrace.internal0 = hd(010) ;
sactrace.t0        = hd(011) ;
sactrace.t1        = hd(012) ;
sactrace.t2        = hd(013) ;
sactrace.t3        = hd(014) ;
sactrace.t4        = hd(015) ;
sactrace.t5        = hd(016) ;
sactrace.t6        = hd(017) ;
sactrace.t7        = hd(018) ;
sactrace.t8        = hd(019) ;
sactrace.t9        = hd(020) ;
sactrace.f         = hd(021) ;
sactrace.resp0     = hd(022) ;
sactrace.resp1     = hd(023) ;
sactrace.resp2     = hd(024) ;
sactrace.resp3     = hd(025) ;
sactrace.resp4     = hd(026) ;
sactrace.resp5     = hd(027) ;
sactrace.resp6     = hd(028) ;
sactrace.resp7     = hd(029) ;
sactrace.resp8     = hd(030) ;
sactrace.resp9     = hd(031) ;
sactrace.stla      = hd(032) ;
sactrace.stlo      = hd(033) ;
sactrace.stel      = hd(034) ;
sactrace.stdp      = hd(035) ;
sactrace.evla      = hd(036) ;
sactrace.evlo      = hd(037) ;
sactrace.evel      = hd(038) ;
sactrace.evdp      = hd(039) ;
sactrace.mag       = hd(040) ;
sactrace.user0     = hd(041) ;
sactrace.user1     = hd(042) ;
sactrace.user2     = hd(043) ;
sactrace.user3     = hd(044) ;
sactrace.user4     = hd(045) ;
sactrace.user5     = hd(046) ;
sactrace.user6     = hd(047) ;
sactrace.user7     = hd(048) ;
sactrace.user8     = hd(049) ;
sactrace.user9     = hd(050) ;
sactrace.dist      = hd(051) ;
sactrace.az        = hd(052) ;
sactrace.baz       = hd(053) ;
sactrace.gcarc     = hd(054) ;
sactrace.internal1 = hd(055) ;
sactrace.internal2 = hd(056) ;
sactrace.depmen    = hd(057) ;
sactrace.cmpaz     = hd(058) ;
sactrace.cmpinc    = hd(059) ;
sactrace.xminimum  = hd(060) ;
sactrace.xmaximum  = hd(061) ;
sactrace.yminimum  = hd(062) ;
sactrace.ymaximum  = hd(063) ;
sactrace.unused1   = hd(064) ;
sactrace.unused2   = hd(065) ;
sactrace.unused3   = hd(066) ;
sactrace.unused4   = hd(067) ;
sactrace.unused5   = hd(068) ;
sactrace.unused6   = hd(069) ;
sactrace.unused7   = hd(070) ;

% convert the integer header

sactrace.nzyear    = ihd(001) ;
sactrace.nzjday    = ihd(002) ;
sactrace.nzhour    = ihd(003) ;
sactrace.nzmin     = ihd(004) ;
sactrace.nzsec     = ihd(005) ;
sactrace.nzmsec    = ihd(006) ;
sactrace.nvhdr     = ihd(007) ;
sactrace.norid     = ihd(008) ;
sactrace.nevid     = ihd(009) ;
sactrace.npts      = ihd(010) ;
sactrace.internal3 = ihd(011) ;
sactrace.nwfid     = ihd(012) ;
sactrace.nxsize    = ihd(013) ;
sactrace.nysize    = ihd(014) ;
sactrace.unused8   = ihd(015) ;
sactrace.iftype    = ihd(016) ;
sactrace.idep      = ihd(017) ;
sactrace.iztype    = ihd(018) ;
sactrace.unused9   = ihd(019) ;
sactrace.iinst     = ihd(020) ;
sactrace.istreg    = ihd(021) ;
sactrace.ievreg    = ihd(022) ;
sactrace.ievtyp    = ihd(023) ;
sactrace.iqual     = ihd(024) ;
sactrace.isynth    = ihd(025) ;
sactrace.imagtyp   = ihd(026) ;
sactrace.imagsrc   = ihd(027) ;
sactrace.unused10  = ihd(028) ;
sactrace.unused11  = ihd(029) ;
sactrace.unused12  = ihd(030) ;
sactrace.unused13  = ihd(031) ;
sactrace.unused14  = ihd(032) ;
sactrace.unused15  = ihd(033) ;
sactrace.unused16  = ihd(034) ;
sactrace.unused17  = ihd(035) ;
sactrace.leven     = ihd(036) ;
sactrace.lpspol    = ihd(037) ;
sactrace.lovrok    = ihd(038) ;
sactrace.lcalda    = ihd(039) ;
sactrace.unused18  = ihd(040) ;

% get the character header
sactrace.kstnm  = [fread(fid,8,'8*char=>char')'] ;
sactrace.kevnm  = [fread(fid,16,'16*char=>char')'] ;
sactrace.khole  = [fread(fid,8,'8*char=>char')'] ;
sactrace.ko     = [fread(fid,8,'8*char=>char')'] ;
sactrace.ka     = [fread(fid,8,'8*char=>char')'] ;
sactrace.kt0    = [fread(fid,8,'8*char=>char')'] ;
sactrace.kt1    = [fread(fid,8,'8*char=>char')'] ;
sactrace.kt2    = [fread(fid,8,'8*char=>char')'] ;
sactrace.kt3    = [fread(fid,8,'8*char=>char')'] ;
sactrace.kt4    = [fread(fid,8,'8*char=>char')'] ;
sactrace.kt5    = [fread(fid,8,'8*char=>char')'] ;
sactrace.kt6    = [fread(fid,8,'8*char=>char')'] ;
sactrace.kt7    = [fread(fid,8,'8*char=>char')'] ;
sactrace.kt8    = [fread(fid,8,'8*char=>char')'] ;
sactrace.kt9    = [fread(fid,8,'8*char=>char')'] ;
sactrace.kf     = [fread(fid,8,'8*char=>char')'] ;
sactrace.kuser0 = [fread(fid,8,'8*char=>char')'] ;
sactrace.kuser1 = [fread(fid,8,'8*char=>char')'] ;
sactrace.kuser2 = [fread(fid,8,'8*char=>char')'] ;
sactrace.kcmpnm = [fread(fid,8,'8*char=>char')'] ;
sactrace.knetwk = [fread(fid,8,'8*char=>char')'] ;
sactrace.kdatrd = [fread(fid,8,'8*char=>char')'] ;
sactrace.kinst  = [fread(fid,8,'8*char=>char')'] ;

sactrace.endian = tmpEndian;

% read the trace
sactrace.x1 = fread(fid,sactrace.npts,'float32') ;

if (length(find(sactrace.iftype==[2 3 4 51]))==1) 
   sactrace.x2 = fread(fid,sactrace.npts,'float32') ;
end
if (sactrace.iftype==51) 
   sactrace.x3 = fread(fid,sactrace.npts,'float32') ;
end

fclose(fid) ;

% end of MSAC_READ.M

