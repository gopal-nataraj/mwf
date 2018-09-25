  function [sp_grad, sm_grad] = dess_2comp_gradx(M0, ff, T1f, T1s, T2f, T2s,...
      kap, Dwf, Dws, flip, TR, TEp, TEm, varargin)
%|function [sp_grad, sm_grad] = dess_2comp_gradx(M0, ff, T1f, T1s, T2f, T2s,...
%|    kap, Dwf, Dws, flip, TR, TEp, TEm, varargin)
%|
%|  2-compartment dess signal gradients
%|    without exchange
%|      uses matlab-generated subfunctions 
%|    with exchange
%|      constructs gradient with calls to dess_2comp_deriv(...)
%|      neglects compartmental exchange up until/prior to echo times
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
%|    kfs       [(odims)]       fast -> slow exchange rate (kHz)                def: zeros(odims)
%|    ksf       [(odims)]       slow -> fast exchange rate (kHz)                def: kfs*(ff/fs)
%|    R2p_f     [(odims)]       fast component broadening linewidth (kHz)       def: zeros(odims)
%|    R2p_s     [(odims)]       slow component broadening linewidth (kHz)       def: zeros(odims)
%|    exchg     false|true      model exchange effects                          def: false
%|    mag       false|true      toggle magnitude signal off|on                  def: true
%|
%|  outputs
%|    sp_grad   [(odims) L]     dess 'defocusing' echo signal gradient at time t=TEp after RF
%|    sm_grad   [(odims) L]     dess 'refocusing' echo signal gradient at time t=TEm before RF
%|  
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2016-04-15      original
%|    1.2       2016-08-15      added magnitude signal option
%|    1.3       2016-08-17      exchg only turned on if nonzero kfs specified
%|    2.1       2018-01-05      now requires separate subroutines for magnitude signal case
%|    2.2       2018-01-19      now requires separate subroutines if comp off-res frequencies are equal
%|    2.3       2018-01-21      now omits fs as an independent variable (forcing fs := 1-ff)

% default values
arg.mask = [];
arg.kfs = [];
arg.ksf = [];
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
  N = prod(odims);
else
  N = numel(arg.mask(arg.mask));
end

% if kfs specified, make sure exchg is on
% if nonzero kfs specified, make sure exchg is off and set to zero
if isempty(arg.kfs)
  if arg.exchg
    warn('Exchange modeled but not specified? kfs/ksf set to zero.');
  end
  arg.kfs = zeros(odims);
elseif ~arg.exchg && nnz(arg.kfs(arg.mask))>0
  warn('Since nonzero kfs set, now modeling with exchange.');
  arg.exchg = true;
end

% if no ksf specified, assume two-compartment equilibrium
if isempty(arg.ksf)
  tmp = bsxfun(@minus, 1, ff);
  arg.ksf = arg.kfs .* div0(ff,tmp);
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
kfs = masker(arg.kfs, arg.mask);
ksf = masker(arg.ksf, arg.mask);
kap = masker(kap, arg.mask);
Dwf = masker(Dwf, arg.mask);
Dws = masker(Dws, arg.mask);
R2pf = masker(arg.R2p_f, arg.mask);
R2ps = masker(arg.R2p_s, arg.mask);

% including exchange, gradients must be constructed numerically
if arg.exchg  
  spxy_tep_rgrad = NaN(N, 7);
  smxy_tem_rgrad = NaN(N, 7);
  [spxy_tep_rgrad(:,1), smxy_tem_rgrad(:,1)] = dess_2comp_deriv(...
    M0, ff, T1f, T1s, T2f, T2s, kfs, kap, Dwf, Dws, flip, TR, TEp, TEm, 'M0',...
    'fs', fs, 'ksf', ksf, 'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
  [spxy_tep_rgrad(:,2), smxy_tem_rgrad(:,2)] = dess_2comp_deriv(...
    M0, ff, T1f, T1s, T2f, T2s, kfs, kap, Dwf, Dws, flip, TR, TEp, TEm, 'ff',...
    'fs', fs, 'ksf', ksf, 'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
  [spxy_tep_rgrad(:,3), smxy_tem_rgrad(:,3)] = dess_2comp_deriv(...
    M0, ff, T1f, T1s, T2f, T2s, kfs, kap, Dwf, Dws, flip, TR, TEp, TEm, 'T1f',...
    'fs', fs, 'ksf', ksf, 'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
  [spxy_tep_rgrad(:,4), smxy_tem_rgrad(:,4)] = dess_2comp_deriv(...
    M0, ff, T1f, T1s, T2f, T2s, kfs, kap, Dwf, Dws, flip, TR, TEp, TEm, 'T1s',...
    'fs', fs, 'ksf', ksf, 'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
  [spxy_tep_rgrad(:,5), smxy_tem_rgrad(:,5)] = dess_2comp_deriv(...
    M0, ff, T1f, T1s, T2f, T2s, kfs, kap, Dwf, Dws, flip, TR, TEp, TEm, 'T2f',...
    'fs', fs, 'ksf', ksf, 'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
  [spxy_tep_rgrad(:,6), smxy_tem_rgrad(:,6)] = dess_2comp_deriv(...
    M0, ff, T1f, T1s, T2f, T2s, kfs, kap, Dwf, Dws, flip, TR, TEp, TEm, 'T2s',...
    'fs', fs, 'ksf', ksf, 'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
  [spxy_tep_rgrad(:,7), smxy_tem_rgrad(:,7)] = dess_2comp_deriv(...
    M0, ff, T1f, T1s, T2f, T2s, kfs, kap, Dwf, Dws, flip, TR, TEp, TEm, 'kfs',...
    'fs', fs, 'ksf', ksf, 'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
% neglecting exchange, use analytical gradients
elseif arg.mag && sum(abs(Dwf-Dws))<eps
  spxy_tep_rgrad = dess_2comp_echo1_exchg0_mag1_freq0_gradx(...
    M0, ff, T1f, T1s, T2f, T2s, kap,...
    flip, TR, TEp, R2pf, R2ps);                         % [N 6] 
  smxy_tem_rgrad = dess_2comp_echo2_exchg0_mag1_freq0_gradx(...
    M0, ff, T1f, T1s, T2f, T2s, kap,...
    flip, TR, TEm, R2pf, R2ps);                         % [N 6]
elseif arg.mag
  error('todo: gradients of magnitude signals with diff comp off-res frequencies');
else
  error('todo: gradients of complex signals');
end

% protect against unsafe division (0/0 set to 0)
spxy_tep_rgrad(isnan(spxy_tep_rgrad)) = 0;
smxy_tem_rgrad(isnan(smxy_tem_rgrad)) = 0;

% set infinities to zero for safety
spxy_tep_rgrad(isinf(spxy_tep_rgrad)) = 0;
smxy_tem_rgrad(isinf(smxy_tem_rgrad)) = 0;

% embed back into original array size
sp_grad = embed(spxy_tep_rgrad, arg.mask);              % [(odims) L]  
sm_grad = embed(smxy_tem_rgrad, arg.mask);              % [(odims) L]  
end
