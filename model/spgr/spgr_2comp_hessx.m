  function [ss_hess] = spgr_2comp_hessx(M0, ff, T1f, T1s, T2f, T2s,...
    kap, Dwf, Dws, flip, TR, TE, varargin)
%|function [ss_hess] = spgr_2comp_hessx(M0, ff, T1f, T1s, T2f, T2s,...
%|  kap, Dwf, Dws, flip, TR, TE, varargin)
%|
%|  2-compartment spgr signal hessian
%|    without exchange
%|      computes with exact analytical expression
%|    with exchange
%|      todo
%|
%|  inputs
%|    M0        [(odims)]       spin density
%|    ff        [(odims)]       fast component spin fraction
%|    T1f       [(odims)]       fast component spin-lattice relaxation time     ms
%|    T1s       [(odims)]       slow component spin-lattice relaxation time     ms
%|    T2f       [(odims)]       fast component spin-spin relaxation time        ms
%|    T2s       [(odims)]       slow component spin-spin relaxation time        ms
%|    kap       [(odims)]       flip angle scaling                  
%|    Dwf       [(odims)]       fast component off-resonance field              kHz
%|    Dws       [(odims)]       slow component off-resonance field              kHz
%|    flip      [1]             nominal nutation angle                          rad
%|    TR        [1]             repetition time                                 ms
%|    TE        [1]             echo time                                       ms
%|    
%|  options
%|    mask      [(odims)]       object mask                                     def: true(odims)
%|    R2p_f     [(odims)]       fast component broadening linewidth (kHz)       def: zeros(odims)
%|    R2p_s     [(odims)]       slow component broadening linewidth (kHz)       def: zeros(odims)
%|    exchg     false|true      model exchange effects                          def: false
%|    mag       false|true      toggle magnitude signal off|on                  def: true
%|
%|  outputs
%|    ss_hess   [(odims) L L]   spgr signal hessian at time t=TE after RF
%|
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2016-10-26      original
%|    2.1       2018-01-21      now requires separate subroutines for magnitude signal case
%|                              now requires separate subroutines if comp off-res frequencies are equal
%|                              now omits fs as an independent variable (forcing fs := 1-ff)

% default values
arg.mask = [];
arg.R2p_f = [];
arg.R2p_s = [];
arg.exchg = false;
arg.mag = true;

% substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% dimensions
odims = size(ff);
 
% if no mask specified, extrapolate to all voxels
if isempty(arg.mask)
  arg.mask = true(odims);
end

% if no fast-component broadening linewidth specified, set to zero
if isempty(arg.R2p_f)
  arg.R2p_f = zeros(odims);
end

% if no slow-component broadening linewidth specified, set to zero
if isempty(arg.R2p_s)
  arg.R2p_s = zeros(odims);
end

% vectorize inputs
M0 = masker(M0, arg.mask);
ff = masker(ff, arg.mask);
T1f = masker(T1f, arg.mask);
T1s = masker(T1s, arg.mask);
T2f = masker(T2f, arg.mask);
T2s = masker(T2s, arg.mask);
kap = masker(kap, arg.mask);
Dwf = masker(Dwf, arg.mask);
Dws = masker(Dws, arg.mask);
R2pf = masker(arg.R2p_f, arg.mask);
R2ps = masker(arg.R2p_s, arg.mask);

if arg.exchg
  error('todo: hessian of signal with physical exchange');
elseif arg.mag && sum(abs(Dwf-Dws))<eps
  ssxy_te_hess = spgr_2comp_exchg0_mag1_freq0_hessx(...
    M0, ff, T1f, T1s, T2f, T2s,...
    kap, flip, TR, TE, R2pf, R2ps);           % [N L L]
elseif arg.mag
  error('todo: hessian of magnitude signal with diff comp off-res frequencies');
else
  error('todo: hessian of complex signal');
end

% protect against unsafe division (0/0 set to 0)
ssxy_te_hess(isnan(ssxy_te_hess)) = 0;

% set infinities to zero for safety
ssxy_te_hess(isinf(ssxy_te_hess)) = 0;

% embed back into original array size
ss_hess = embed(ssxy_te_hess, arg.mask);      % [(odims) L L]
end
  