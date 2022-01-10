% MSAC toolkit function.
%
%  [] = msac_write(filename,sac_trace_structure) ;
%
%  Write a (single) SAC trace structure to a file. See msac_new.m for
%  details of how to create a SAC structure. The value in the field
%  sac_trace_structure.endian is used to determine the byte-order of
%  the resulting file, if it exists, otherwise big-endian is used by 
%  default

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

function msac_write(str,sactrace) ;

if (isfield(sactrace,'endian'))
   [fid,message] = fopen(str,'wb',sactrace.endian) ;
else
   [fid,message] = fopen(str,'wb','b') ;
end
   
if (~isempty(message))
   fprintf('Cannot open file for writing: ''%s''\n',str) ;
   return ;
end

hd(001)  = sactrace.delta    ;
hd(002)  = sactrace.depmin   ;
hd(003)  = sactrace.depmax   ;
hd(004)  = sactrace.scale    ;
hd(005)  = sactrace.odelta   ;
hd(006)  = sactrace.b        ;
hd(007)  = sactrace.e        ;
hd(008)  = sactrace.o        ;
hd(009)  = sactrace.a        ;
hd(010)  = sactrace.internal0;
hd(011)  = sactrace.t0       ;
hd(012)  = sactrace.t1       ;
hd(013)  = sactrace.t2       ;
hd(014)  = sactrace.t3       ;
hd(015)  = sactrace.t4       ;
hd(016)  = sactrace.t5       ;
hd(017)  = sactrace.t6       ;
hd(018)  = sactrace.t7       ;
hd(019)  = sactrace.t8       ;
hd(020)  = sactrace.t9       ;
hd(021)  = sactrace.f        ;
hd(022)  = sactrace.resp0    ;
hd(023)  = sactrace.resp1    ;
hd(024)  = sactrace.resp2    ;
hd(025)  = sactrace.resp3    ;
hd(026)  = sactrace.resp4    ;
hd(027)  = sactrace.resp5    ;
hd(028)  = sactrace.resp6    ;
hd(029)  = sactrace.resp7    ;
hd(030)  = sactrace.resp8    ;
hd(031)  = sactrace.resp9    ;
hd(032)  = sactrace.stla     ;
hd(033)  = sactrace.stlo     ;
hd(034)  = sactrace.stel     ;
hd(035)  = sactrace.stdp     ;
hd(036)  = sactrace.evla     ;
hd(037)  = sactrace.evlo     ;
hd(038)  = sactrace.evel     ;
hd(039)  = sactrace.evdp     ;
hd(040)  = sactrace.mag      ;
hd(041)  = sactrace.user0    ;
hd(042)  = sactrace.user1    ;
hd(043)  = sactrace.user2    ;
hd(044)  = sactrace.user3    ;
hd(045)  = sactrace.user4    ;
hd(046)  = sactrace.user5    ;
hd(047)  = sactrace.user6    ;
hd(048)  = sactrace.user7    ;
hd(049)  = sactrace.user8    ;
hd(050)  = sactrace.user9    ;
hd(051)  = sactrace.dist     ;
hd(052)  = sactrace.az       ;
hd(053)  = sactrace.baz      ;
hd(054)  = sactrace.gcarc    ;
hd(055)  = sactrace.internal1;
hd(056)  = sactrace.internal2;
hd(057)  = sactrace.depmen   ;
hd(058)  = sactrace.cmpaz    ;
hd(059)  = sactrace.cmpinc   ;
hd(060)  = sactrace.xminimum ;
hd(061)  = sactrace.xmaximum ;
hd(062)  = sactrace.yminimum ;
hd(063)  = sactrace.ymaximum ;
hd(064)  = sactrace.unused1  ;
hd(065)  = sactrace.unused2  ;
hd(066)  = sactrace.unused3  ;
hd(067)  = sactrace.unused4  ;
hd(068)  = sactrace.unused5  ;
hd(069)  = sactrace.unused6  ;
hd(070)  = sactrace.unused7  ;

ihd(001) = sactrace.nzyear    ;
ihd(002) = sactrace.nzjday    ;
ihd(003) = sactrace.nzhour    ;
ihd(004) = sactrace.nzmin     ;
ihd(005) = sactrace.nzsec     ;
ihd(006) = sactrace.nzmsec    ;
ihd(007) = sactrace.nvhdr     ;
ihd(008) = sactrace.norid     ;
ihd(009) = sactrace.nevid     ;
ihd(010) = sactrace.npts      ;
ihd(011) = sactrace.internal3 ;
ihd(012) = sactrace.nwfid     ;
ihd(013) = sactrace.nxsize    ;
ihd(014) = sactrace.nysize    ;
ihd(015) = sactrace.unused8   ;
ihd(016) = sactrace.iftype    ;
ihd(017) = sactrace.idep      ;
ihd(018) = sactrace.iztype    ;
ihd(019) = sactrace.unused9   ;
ihd(020) = sactrace.iinst     ;
ihd(021) = sactrace.istreg    ;
ihd(022) = sactrace.ievreg    ;
ihd(023) = sactrace.ievtyp    ;
ihd(024) = sactrace.iqual     ;
ihd(025) = sactrace.isynth    ;
ihd(026) = sactrace.imagtyp   ;
ihd(027) = sactrace.imagsrc   ;
ihd(028) = sactrace.unused10  ;
ihd(029) = sactrace.unused11  ;
ihd(030) = sactrace.unused12  ;
ihd(031) = sactrace.unused13  ;
ihd(032) = sactrace.unused14  ;
ihd(033) = sactrace.unused15  ;
ihd(034) = sactrace.unused16  ;
ihd(035) = sactrace.unused17  ;
ihd(036) = sactrace.leven     ;
ihd(037) = sactrace.lpspol    ;
ihd(038) = sactrace.lovrok    ;
ihd(039) = sactrace.lcalda    ;
ihd(040) = sactrace.unused18  ;


% output the float header
fwrite(fid,hd,'float32') ;

% output the integer header
fwrite(fid,ihd,'int') ;

nullstr8 =  char([0 0 0 0 0 0 0 0]) ;
nullstr16 = char([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]) ;

% output the character header, padding with nulls where necessary
out=[sactrace.kstnm  nullstr8]; fwrite(fid,out(1:8), 'char') ;
out=[sactrace.kevnm  nullstr16];fwrite(fid,out(1:16),'char') ;
out=[sactrace.khole  nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.ko     nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.ka     nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kt0    nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kt1    nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kt2    nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kt3    nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kt4    nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kt5    nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kt6    nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kt7    nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kt8    nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kt9    nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kf     nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kuser0 nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kuser1 nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kuser2 nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kcmpnm nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.knetwk nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kdatrd nullstr8]; fwrite(fid,out(1:8),'char' ) ;
out=[sactrace.kinst  nullstr8]; fwrite(fid,out(1:8),'char' ) ;

% write the trace
fwrite(fid,sactrace.x1,'float32') ;

if (length(find(sactrace.iftype==[2 3 4 51]))==1) 
   fwrite(fid,sactrace.x2,'float32') ;
end
if (sactrace.iftype==51) 
   fwrite(fid,sactrace.x3,'float32') ;
end

fclose(fid) ;

% end of MSAC_WRITE.M

