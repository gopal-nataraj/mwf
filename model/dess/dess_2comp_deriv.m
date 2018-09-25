  function [sp_der, sm_der] = dess_2comp_deriv(M0, ff, T1f, T1s, T2f, T2s, kfs,...
      kap, Dwf, Dws, flip, TR, TEp, TEm, ind_var, varargin)
%|function [sp_der, sm_der] = dess_2comp_deriv(M0, ff, T1f, T1s, T2f, T2s, kfs,...
%|    kap, Dwf, Dws, flip, TR, TEp, TEm, ind_var, varargin)
%|
%|  2-compartment dess signal numerical partial derivatives
%|    neglects compartmental exchange up until/prior to echo times
%|    warning: derivative accuracy vs. numerical stability depends strongly on delta!
%|
%|  inputs
%|    M0        [(odims)]       spin density
%|    ff        [(odims)]       fast component spin fraction
%|    T1f       [(odims)]       fast component spin-lattice relaxation time     ms
%|    T1s       [(odims)]       slow component spin-lattice relaxation time     ms
%|    T2f       [(odims)]       fast component spin-spin relaxation time        ms
%|    T2s       [(odims)]       slow component spin-spin relaxation time        ms
%|    kfs       [(odims)]       fast -> slow exchange rate                      kHz
%|    kap       [(odims)]       flip angle scaling                  
%|    Dwf       [(odims)]       fast component off-resonance field              kHz
%|    Dws       [(odims)]       slow component off-resonance field              kHz
%|    flip      [1]             nominal nutation angle                          rad
%|    TR        [1]             repetition time                                 ms
%|    TEp       [1]             'defocusing' echo time                          ms
%|    TEm       [1]             'refocusing' echo time                          ms
%|    ind_var   '1'             variable by which to take partial derivative
%|
%|  options
%|    mask      [(odims)]       object mask                                     def: true(odims)
%|    fs        [(odims)]       slow component spin fraction                    def: 1-ff
%|    ksf       [(odims)]       slow -> fast exchange rate (kHz)                def: kfs*(ff/fs)
%|    R2p_f     [(odims)]       fast component broadening linewidth (kHz)       def: zeros(odims)
%|    R2p_s     [(odims)]       slow component broadening linewidth (kHz)       def: zeros(odims)
%|    delta     [1]             proportion by which to perturb signal           def: 100*eps
%|    exchg     false|true      model exchange off|on                           def: false
%|  
%|  outputs
%|    sp_der    [(odims)]       'defocusing' echo numerical partial derivative
%|    sm_der    [(odims)]       'refocusing' echo numerical partial derivative
%|  
%|  copyright 2016, gopal nataraj, university of michigan 
%|        
%|  version control  
%|    1.1       2016-03-25      original
%|    1.2       2016-03-28      changed from 9->7 independent variables, per constraint
%|    1.3       2016-03-30      added optional arguments, flip angle spatial variations
%|    1.4       2016-04-15      modified calls to dess_2comp(...), reduced delta greatly

% default values
arg.mask = [];
arg.fs = [];
arg.ksf = [];
arg.R2p_f = [];
arg.R2p_s = [];
arg.delta = 100*eps;
arg.exchg = false;

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

% if no fs specified, assume two compartments only
if isempty(arg.fs)
  arg.fs = 1 - ff;
end
ff = min(ff, 1-eps);
arg.fs = max(arg.fs, eps);

% if no ksf specified, assume two-compartment equilibrium
if isempty(arg.ksf)
  arg.ksf = kfs .* div0(ff,arg.fs);
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
fs = masker(arg.fs, arg.mask);
T1f = masker(T1f, arg.mask);
T1s = masker(T1s, arg.mask);
T2f = masker(T2f, arg.mask);
T2s = masker(T2s, arg.mask);
kfs = masker(kfs, arg.mask);
ksf = masker(arg.ksf, arg.mask);
kap = masker(kap, arg.mask);
Dwf = masker(Dwf, arg.mask);
Dws = masker(Dws, arg.mask);
R2pf = masker(arg.R2p_f, arg.mask);
R2ps = masker(arg.R2p_s, arg.mask);

% compute perturbed signal
switch ind_var
case 'M0'
  dvar = arg.delta * M0;
  [sp2, sm2] = dess_2comp(M0+dvar, ff, T1f, T1s, T2f, T2s,...
    kap, Dwf, Dws, flip, TR, TEp, TEm,...
    'fs', fs, 'kfs', kfs, 'ksf', ksf,...
    'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
case 'ff'
  dvar = arg.delta * ff;
  [sp2, sm2] = dess_2comp(M0, ff+dvar, T1f, T1s, T2f, T2s,...
    kap, Dwf, Dws, flip, TR, TEp, TEm,...
    'fs', fs, 'kfs', kfs, 'ksf', ksf,...
    'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
case 'T1f'
  dvar = arg.delta * T1f;
  [sp2, sm2] = dess_2comp(M0, ff, T1f+dvar, T1s, T2f, T2s,...
    kap, Dwf, Dws, flip, TR, TEp, TEm,...
    'fs', fs, 'kfs', kfs, 'ksf', ksf,...
    'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
case 'T1s'
  dvar = arg.delta * T1s;
  [sp2, sm2] = dess_2comp(M0, ff, T1f, T1s+dvar, T2f, T2s,...
    kap, Dwf, Dws, flip, TR, TEp, TEm,...
    'fs', fs, 'kfs', kfs, 'ksf', ksf,...
    'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
case 'T2f'
  dvar = arg.delta * T2f;
  [sp2, sm2] = dess_2comp(M0, ff, T1f, T1s, T2f+dvar, T2s,...
    kap, Dwf, Dws, flip, TR, TEp, TEm,...
    'fs', fs, 'kfs', kfs, 'ksf', ksf,...
    'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
case 'T2s'
  dvar = arg.delta * T2s;
  [sp2, sm2] = dess_2comp(M0, ff, T1f, T1s, T2f, T2s+dvar,...
    kap, Dwf, Dws, flip, TR, TEp, TEm,...
    'fs', fs, 'kfs', kfs, 'ksf', ksf,...
    'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
case 'kfs'
  dvar = arg.delta * kfs;
  [sp2, sm2] = dess_2comp(M0, ff, T1f, T1s, T2f, T2s,...
    kap, Dwf, Dws, flip, TR, TEp, TEm,...
    'fs', fs, 'kfs', kfs+dvar, 'ksf', ksf,...
    'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
case 'kap'
  dvar = arg.delta * kap;
  [sp2, sm2] = dess_2comp(M0, ff, T1f, T1s, T2f, T2s,...
    kap+dvar, Dwf, Dws, flip, TR, TEp, TEm,...
    'fs', fs, 'kfs', kfs, 'ksf', ksf,...
    'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
case 'Dwf'
  dvar = arg.delta * Dwf;
  [sp2, sm2] = dess_2comp(M0, ff, T1f, T1s, T2f, T2s,...
    kap, Dwf+dvar, Dws, flip, TR, TEp, TEm,...
    'fs', fs, 'kfs', kfs, 'ksf', ksf,...
    'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
case 'Dws'
  dvar = arg.delta * Dws;
  [sp2, sm2] = dess_2comp(M0, ff, T1f, T1s, T2f, T2s,...
    kap, Dwf, Dws+dvar, flip, TR, TEp, TEm,...
    'fs', fs, 'kfs', kfs, 'ksf', ksf,...
    'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
otherwise
  error('derivative failed: unknown variable ''%s''.', ind_var);
end

% compute baseline signal
[sp1, sm1] = dess_2comp(M0, ff, T1f, T1s, T2f, T2s,...
  kap, Dwf, Dws, flip, TR, TEp, TEm,...
  'fs', fs, 'kfs', kfs, 'ksf', ksf,... 
  'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
        
% numerical derivatives
% embed back into original array sizes
sp_der = embed(div0(sp2-sp1, dvar), arg.mask);
sm_der = embed(div0(sm2-sm1, dvar), arg.mask);
  end