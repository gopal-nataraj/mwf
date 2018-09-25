  function [sp_gx_gp, sm_gx_gp] = dess_2comp_gradx_gradp(M0, ff, T1f, T1s, T2f, T2s,...
      kap, Dwf, Dws, flip, TR, TEp, TEm, varargin)
%|function [sp_gx_gp, sm_gx_gp] = dess_2comp_gradx_gradp(M0, ff, T1f, T1s, T2f, T2s,...
%|    kap, Dwf, Dws, flip, TR, TEp, TEm, varargin)
%|
%|  2-compartment dess mixed partial derivatives
%|    without exchange
%|      uses matlab-generated subfunctions dess_2comp_echo*_exchg0_gradx_gradp_gen(...)
%|    todo: with exchange
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
%|    TEp       [1]             'defocusing' echo time                          ms
%|    TEm       [1]             'refocusing' echo time                          ms
%|    
%|  options
%|    mask      [(odims)]       object mask                                     def: true(odims)
%|    R2p_f     [(odims)]       fast component broadening linewidth (kHz)       def: zeros(odims)
%|    R2p_s     [(odims)]       slow component broadening linewidth (kHz)       def: zeros(odims)
%|    exchg     false|true      model exchange effects                          def: false
%|    mag       false|true      toggle magnitude signal off|on                  def: true
%|
%|  outputs
%|    sp_gx_gp  [(odims) L P]   dess 'defocusing' echo mixed second partial derivatives
%|    sm_gx_gp  [(odims) L P]   dess 'refocusing' echo mixed second partial derivatives
%| 
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2016-04-21      original
%|    2.1       2018-01-05      now requires separate subroutines for magnitude signal case
%|    2.2       2018-01-19      now requires separate subroutines if comp off-res frequencies are equal
%|    2.3       2018-01-21      now omits fs as an independent variable (forcing fs := 1-ff)

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
  error('todo: mixed analytical gradients with exchange.');
elseif arg.mag && sum(abs(Dwf-Dws))<eps
  sp_gx_gp = dess_2comp_echo1_exchg0_mag1_freq0_gradx_gradp(...
    M0, ff, T1f, T1s, T2f, T2s,...
    kap, flip, TR, TEp, R2pf, R2ps);                    % [N L P]
  sm_gx_gp = dess_2comp_echo2_exchg0_mag1_freq0_gradx_gradp(...
    M0, ff, T1f, T1s, T2f, T2s,...
    kap, flip, TR, TEm, R2pf, R2ps);                    % [N L P]
elseif arg.mag
  error('todo: mixed gradients of magnitude signals with diff comp off-res frequencies');
else
  error('todo: mixed gradients of complex signals');
end

% protect against unsafe division (0/0 set to 0)
sp_gx_gp(isnan(sp_gx_gp)) = 0;
sm_gx_gp(isnan(sm_gx_gp)) = 0;

% set infinities to zero for safety
sp_gx_gp(isinf(sp_gx_gp)) = 0;
sm_gx_gp(isinf(sm_gx_gp)) = 0;

% embed back into original array size
sp_gx_gp = embed(sp_gx_gp, arg.mask);                   % [(odims) L P]
sm_gx_gp = embed(sm_gx_gp, arg.mask);                   % [(odims) L P]
end
