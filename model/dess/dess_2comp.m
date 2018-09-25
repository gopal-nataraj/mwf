  function [sp, sm] = dess_2comp(M0, ff, T1f, T1s, T2f, T2s,...
    kap, Dwf, Dws, flip, TR, TEp, TEm, varargin)
%|function [sp, sm] = dess_2comp(M0, ff, T1f, T1s, T2f, T2s,...
%|  kap, Dwf, Dws, flip, TR, TEp, TEm, varargin)
%|
%|  2-compartment noiseless dess signal models, with optional first-order exchange 
%|  without exchange
%|    computes with exact analytical expressions
%|  with exchange
%|    requires numerical integration over off-resonance phases
%|    neglects compartmental exchange up until/prior to echo times
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
%|    ncyc      [1]             number of unbalanced gradient spoiling cycles   def: 1
%|                              if not an integer, models partial spoiling
%|    nph       [1]             number of off-resonant phases to sample         def: 100
%|    exchg     false|true      toggle exchange off|on                          def: false
%|    mag       false|true      toggle magnitude signal off|on                  def: true
%|    
%|  outputs  
%|    sp        [(odims)]       'defocused' dess signal at t=TEp after RF
%|    sm        [(odims)]       'refocused' dess signal at t=TEm before RF
%|    
%|  copyright 2016, gopal nataraj, university of michigan 
%|        
%|  version control  
%|    1.1       2016-03-23      original
%|    1.2       2016-03-28      added default fs, ksf values as per constraints
%|    1.3       2016-03-30      added flip angle spatial variation
%|    1.4       2016-04-15      added analytical options for without exchange
%|    1.5       2016-08-15      added magnitude signal option
%|    1.6       2016-08-17      exchg only turned on if nonzero kfs specified
%|    2.1       2018-01-21      now omits fs as an independent variable (forcing fs := 1-ff)

% default values
arg.mask = [];
arg.kfs = [];
arg.ksf = [];
arg.R2p_f = [];
arg.R2p_s = [];
arg.ncyc = 1;
arg.nph = 100;
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

% with exchange, use numerical integration
if arg.exchg
  % sample spoiling cycles
  if rem(arg.ncyc,1)~=0
    warn('Non-integral number of spoiling cycles.');
  end
  phi = linspace(0, 2*pi*arg.ncyc, arg.nph);

  % compute ss magnetizations at t=0\pm with fixed off-resonant phase
  mp_0 = NaN(6, N, arg.nph);                          % [6 N arg.nph]
  mm_0 = NaN(6, N, arg.nph);                          % [6 N arg.nph]
  for p = 1:arg.nph
    [mp_0(:,:,p), mm_0(:,:,p)] = ss_2comp_ph_fixed(...
      M0, ff, T1f, T1s, T2f, T2s, kfs, ksf, kap, flip, TR, phi(p));
  end

  % numerically integrate/normalize over off-resonant phases
  sp_0 = div0(trapz(phi, mp_0, 3), 2*pi*arg.ncyc);    % [6 N]
  sm_0 = div0(trapz(phi, mm_0, 3), 2*pi*arg.ncyc);    % [6 N]

  % extract complex, compartment-wise signals
  spxy_0_f = col(complex(sp_0(1,:), sp_0(3,:)));      % [N]
  spxy_0_s = col(complex(sp_0(2,:), sp_0(4,:)));      % [N]
  smxy_0_f = col(complex(sm_0(1,:), sm_0(3,:)));      % [N]
  smxy_0_s = col(complex(sm_0(2,:), sm_0(4,:)));      % [N]

  % apply compartment-wise dephasing between tip and signal collection
  % warn: ignores compartmental exchange
  spxy_tep_f = spxy_0_f .* exp(-TEp * (R2pf + 1./T2f)) .* exp(+1i * TEp * Dwf);
  spxy_tep_s = spxy_0_s .* exp(-TEp * (R2ps + 1./T2s)) .* exp(+1i * TEp * Dws);
  smxy_tem_f = smxy_0_f .* exp(-TEm * (R2pf - 1./T2f)) .* exp(-1i * TEm * Dwf);
  smxy_tem_s = smxy_0_s .* exp(-TEm * (R2ps - 1./T2s)) .* exp(-1i * TEm * Dws);

  % combine compartmental contributions
  % embed back into original array sizes
  sp = spxy_tep_f + spxy_tep_s;
  sm = smxy_tem_f + smxy_tem_s;
% without exchange, use analytical models
elseif ~arg.mag
  sp = dess_2comp_echo1_exchg0_mag0_freq0(...
    M0, ff, T1f, T1s, T2f, T2s, kap, Dwf, Dws,...
    flip, TR, TEp, R2pf, R2ps);
  sm = dess_2comp_echo2_exchg0_mag0_freq0(...
    M0, ff, T1f, T1s, T2f, T2s, kap, Dwf, Dws,...
    flip, TR, TEm, R2pf, R2ps);
elseif arg.mag && sum(abs(Dwf-Dws))<eps
  sp = dess_2comp_echo1_exchg0_mag1_freq0(...
    M0, ff, T1f, T1s, T2f, T2s, kap,...
    flip, TR, TEp, R2pf, R2ps);
  sm = dess_2comp_echo2_exchg0_mag1_freq0(...
    M0, ff, T1f, T1s, T2f, T2s, kap,...
    flip, TR, TEm, R2pf, R2ps);
elseif arg.mag
  error('todo: magnitude signals with diff comp off-res frequencies');
end

% protect against unsafe division (0/0 set to 0)
sp(isnan(sp)) = 0;
sm(isnan(sm)) = 0;

% optional: return magnitude signal
if arg.mag
  sp = abs(sp);
  sm = abs(sm);
end

% embed back into original array size
sp = embed(sp, arg.mask);
sm = embed(sm, arg.mask);
end
  
  
  function [mp_0, mm_0] = ss_2comp_ph_fixed(M0,ff,T1f,T1s,T2f,T2s,kfs,ksf,kap,flip,TR,phi)
%|function [mp_0, mm_0] = ss_2comp_ph_fixed(M0,ff,T1f,T1s,T2f,T2s,kfs,ksf,kap,flip,TR,phi)
%|
%|  2-compartment ss magnetization at t=0\pm with fixed off-resonant phase

% dimensions
N = length(M0);

% add spatial variation to nominal flip angle
a = kap * flip;                             % [N 1]

% slow-relaxing fraction
fs = ones(N,1)-ff;

% within-voxel phase primarily due to gradients, not b0 inhomogeneity
Dwf = div0(phi, TR);
Dws = div0(phi, TR);

% compute ss signal voxel-wise
mp_0 = NaN(6, N);                           % [6 N]
mm_0 = NaN(6, N);                           % [6 N]
for n = 1:N
  % assemble bloch matrix
  Axy  = [-1/T2f(n)-kfs(n),   ksf(n),             Dwf,                0;      
          kfs(n),             -1/T2s(n)-ksf(n),   0,                  Dws;
          -Dwf,               0,                  -1/T2f(n)-kfs(n),   ksf(n);
          0,                  -Dws,               kfs(n),             -1/T2s(n)-ksf(n)];
  Az   = [-1/T1f(n)-kfs(n),   ksf(n);
          kfs(n),            -1/T1s(n)-ksf(n)];
  A    = [Axy, zeros(4,2);
          zeros(2,4), Az];

  % exponential matrix
  expA = expm(A * TR);

  % assemble constant matrix
  c = M0(n) * [zeros(4,1); div0(ff(n),T1f(n)); div0(fs(n),T1s(n))];

  % assemble rotation matrix
  R3   = [1,                  0,                  0;
          0,                  cos(a(n)),       sin(a(n));
          0,                  -sin(a(n)),      cos(a(n))];
  R6 = kron(R3, eye(2));

  % steady-state magnetization expressions
  tmp = (expA - eye(6)) * (A \ c);
  mp_0(:,n) = (eye(6) - R6*expA) \ (R6 * tmp);
  mm_0(:,n) = (eye(6) - expA*R6) \ tmp;
%   mp_0(:,i) = ((eye(6) - R6*expA) \ R6) * (expA - eye(6)) * (A \ c);
%   mm_0(:,i) = ((eye(6) - expA*R6) \ (expA - eye(6))) * (A \ c);
%   mm_0(:,i) = R6 \ mp_0(:,i);
end
end
  