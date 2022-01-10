%MSAC_NEW - Return a new SAC structure
%
%  [sac_trace_structure] = msac_new(data,delta) ;
%
%     Returns a new SAC structure with data x [y, and z], and sampling delta.
%     This sets a default b and e of 0 and (nx-1)*delta respectively. All other
%     headers are set to SAC Null values. data must be an n-by-m matrix, where
%     one of n or m must be 1-3 (1 for time-series; 2 for spectral or x,y, 3
%     for xyz). If size(x)=[2 n], and iftype is not explicitly set, a real-imag. 
%     type file is assumeed.
%
%  [sac_trace_structure] = msac_new(timeseries,delta,HEADER,VALUE,...) ;
%     
%     Further arguments are interpreted as header, value pairs: e.g.,
%     [tr] = msac_new(x,delta,'evdp',600) sets the event depth header
%     to be 600 (kilometres). There (obviously) must be an even number of
%     these. 
%
%  NB. It is assumed that the small dimension of the input matrix data is the
%      number of variables. 
%
%  NB. No 'smart' processing of header values is undertaken (e.g., setting b
%  or e will not set the other; setting event and station location will not
%  set gcarc etc etc). 
%
%  NB. Default byte-order is big-endian, even on a little endian machine. To
%      change, set the endian header to 'l'. This is for compatibility with 
%      MacSAC.

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

function [tr] = msac_new(data,delta,varargin) ;

% check input parameters
if nargin < 2 
   error('Not enough arguments for MSAC_NEW')
end
if length(size(data))>2
   error('Argument data must be a vector or a matrix')
end
if ~isscalar(delta)   
   error('Argument delta must be a scalar')
end

%  figure out dimensions of the input data, transpose if necessary
sz = size(data) ;
nvar = min(sz) ;
npts = max(sz) ;
if (sz(1)>sz(2))
   data = data';
end   
   
if (nvar>3)
   error('Only maximum 3 variables can be represented in a SAC file.')
end   

%  setup SAC null values

   SAC_rnull = -12345.0      ;
   SAC_inull = -12345        ;
   SAC_lnull = -12345        ;
   SAC_cnull = '-12345  '    ;

   sachdrs = ...
          {'delta','depmin','depmax','scale','odelta','b','e','o','a',...
           'internal0','t0','t1','t2','t3','t4','t5','t6','t7','t8','t9',...
           'f','resp0','resp1','resp2','resp3','resp4','resp5','resp6',...
           'resp7','resp8','resp9','stla','stlo','stel','stdp','evla',...
           'evlo','evel','evdp','mag','user0','user1','user2','user3',...
           'user4','user5','user6','user7','user8','user9','dist','az',...
           'baz','gcarc','internal1','internal2','depmen','cmpaz',...
           'cmpinc','xminimum','xmaximum','yminimum','ymaximum',...
           'unused1','unused2','unused3','unused4','unused5','unused6',...
           'unused7','nzyear','nzjday','nzhour','nzmin','nzsec',...
           'nzmsec','nvhdr','norid','nevid','npts','internal3','nwfid',...
           'nxsize','nysize','unused8','iftype','idep','iztype',...
           'unused9','iinst','istreg','ievreg','ievtyp','iqual',...
           'isynth','imagtyp','imagsrc','unused10','unused11',...
           'unused12','unused13','unused14','unused15','unused16',...
           'unused17','leven','lpspol','lovrok','lcalda','unused18',...
           'kstnm','kevnm','khole','ko','ka','kt0','kt1','kt2','kt3',...
           'kt4','kt5','kt6','kt7','kt8','kt9','kf','kuser0','kuser1',...
           'kuser2','kcmpnm','knetwk','kdatrd','kinst'} ;
   
   for i=[1:70]
      eval(sprintf('tr.%s = SAC_rnull ;',sachdrs{i})) ;
   end
      
   tr.delta     = delta ;
   tr.b         = 0.0 ;

   for i=[71:105]
       eval(sprintf('tr.%s = SAC_inull ;',sachdrs{i})) ;
   end

   tr.nvhdr     = 6 ;% default
   tr.npts      = npts ;% number of samples
   
   if nvar==1
      tr.iftype = 1 ; % default (time series file)
      tr.x1 = data ;
   elseif nvar==2
      tr.iftype = 2 ; % spectral file (real-imaginary)
      tr.x1 = data(1,:) ;
      tr.x2 = data(2,:) ;
   elseif nvar==3
      tr.iftype = 51 ; % XYZ file
      tr.x1 = data(1,:) ;
      tr.x2 = data(2,:) ;
      tr.x3 = data(3,:) ;
   end    
   
   tr.idep      = 5 ;% default
   tr.iztype    = 9 ;% default
   tr.ievtyp    = 5 ;% default

   tr.leven     = 1 ;% default
   tr.lpspol    = 0 ;
   tr.lovrok    = 1 ;
   tr.lcalda    = 1 ;
   tr.unused18  = 0 ;

   for i=[111:133]
          eval(sprintf('tr.%s = SAC_cnull ;',sachdrs{i})) ;
   end

%  setup the .endian field (used for writing)
   tr.endian = 'b' ; % default

   
% now process the optional arguments (if there are any)
   if nargin == 2 % no, there aren't
      return ;
   end
   
%  check that we have an even number of arguments
   if mod(length(varargin),2)~=0
      error('Header-Argument pairs must be specified') ;
   end
   
%  now process the list
   for i=1:2:length(varargin)-1

%     check that the header is valid
      tmp = cell(1,133) ;
      [tmp{:}] = deal(varargin{i}) ;
      index1 = strcmp(sachdrs,tmp) ;
      index2 = find(index1==1) ;
      if isempty(index2)
         error('Invalid header specified: %s',varargin{i}) ;
      end

%     check that the value specified is valid, and set it
%     REAL HEADERS
      if index2 >= 1 & index2 <= 70
         if ~isnumeric(varargin{i+1})
            error('Header %s must take a numeric value',varargin{i})
         end
         eval(sprintf('tr.%s = varargin{i+1} ;',varargin{i})) ;
         
%     INTEGER HEADERS
      elseif index2 >= 71 & index2 <= 105
         if ~isnumeric(varargin{i+1})
            error('Header %s must be an integer',varargin{i})
         end
         if mod(varargin{i+1},1)~=0
            error('Header %s must be an integer',varargin{i})
         end
         eval(sprintf('tr.%s = varargin{i+1} ;',varargin{i})) ;

%     LOGICAL HEADERS
      elseif index2 >= 106 & index2 <= 110
         if ~islogical(varargin{i+1})
            error('Header %s must take a logical (0 or 1) value',varargin{i})
         end
         eval(sprintf('tr.%s = varargin{i+1} ;',varargin{i})) ;
        
%     CHARACTER HEADERS
      elseif index2 > 110
         if ~ischar(varargin{i+1})
            error('Header %s must take a string value',varargin{i})
         end
         if index2 == 112
            ncharmax = 16 ;
         else
            ncharmax = 8 ;
         end
         
         if length(varargin{i+1}) > ncharmax
            eval(sprintf('tr.%s = varargin{i+1}(1:ncharmax) ;',varargin{i})) ;
            warning('Header value for %s was truncated',varargin{i}) ;
         else
            eval(sprintf('tr.%s = varargin{i+1} ;',varargin{i})) ;
         end
      end
   end
   
%  finally, check that if the iftype was set, it is consistent with the number 
%  of variables. 

if (nvar==1)
   if (tr.iftype~=1), error('Selected IFTYPE incompatible with data provided'),end
elseif (nvar==2)
   if (length(find(tr.iftype==[2 3 4]))~=1)
      error('Selected IFTYPE incompatible with data provided')
   end
elseif (nvar==3)
   if (tr.iftype~=51), error('Selected IFTYPE incompatible with data provided'),end
end
% end of MSAC_NEW

