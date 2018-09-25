  function F = dess_2comp_crb(ff, T1f, T1s, T2f, T2s,...
      kap, Dwf, Dws, flip, TR, TEp, TEm, sigp, sigm, varargin)
%|function F = dess_2comp_crb(ff, T1f, T1s, T2f, T2s,...
%|    kap, Dwf, Dws, flip, TR, TEp, TEm, sigp, sigm, varargin)
%|
%|  2-compartment dess signal fisher information matrix
%|    assumes fs and ksf are fixed via constraints
%|    assumes Dwf, Dws, kap (flip scaling) all known
%|    if exchg on/off, then L=7/6 total latent parameters
%|
%|  inputs
%|    ff        [(odims)]       fast component spin fraction
%|    T1f       [(odims)]       fast component spin-lattice relaxation time     ms
%|    T1s       [(odims)]       slow component spin-lattice relaxation time     ms
%|    T2f       [(odims)]       fast component spin-spin relaxation time        ms
%|    T2s       [(odims)]       slow component spin-spin relaxation time        ms
%|    kap       [(odims)]       flip angle scaling                  
%|    Dwf       [(odims)]       fast component off-resonance field              kHz
%|    Dws       [(odims)]       slow component off-resonance field              kHz
%|    flip      [nf]            nominal nutation angle                          rad
%|    TR        [nf]            repetition time                                 ms
%|    TEp       [nf]            'defocusing' echo time                          ms
%|    TEm       [nf]            'refocusing' echo time                          ms
%|    sigp      [nf]            'defocusing' echo unit-M0 noise std dev
%|    sigm      [nf]            'refocusing' echo unit-M0 noise std dev
%|
%|  options
%|    mask      [(odims)]       object mask                                     def: true(odims)
%|    M0        [(odims)]       spin density (recommend not changing!)          def: ones(odims)
%|    kfs       [(odims)]       fast -> slow exchange rate (kHz)                def: zeros(odims)
%|    ksf       [(odims)]       slow -> fast exchange rate (kHz)                def: kfs*(ff/fs)
%|    R2p_f     [(odims)]       fast component broadening linewidth (kHz)       def: zeros(odims)
%|    R2p_s     [(odims)]       slow component broadening linewidth (kHz)       def: zeros(odims)
%|    exchg     false|true      model exchange effects                          def: false
%| 
%|  outputs
%|    F         [(odims) L L]   fisher information matrix
%|
%|  copyright 2016, gopal nataraj, university of michigan 
%|        
%|  version control  
%|    1.1       2016-03-24      original
%|    1.2       2016-03-28      corrected constrained derivatives, allowed multiple voxels
%|    1.3       2016-03-30      added optional arguments, flip angle spatial variations
%|    1.4       2016-04-15      modified to work with dess_2comp_gradx(...)
%|    1.5       2018-01-21      now omits fs as an independent variable (forcing fs := 1-ff)

% default values
arg.mask = [];  
arg.M0 = [];
arg.kfs = [];
arg.ksf = [];
arg.R2p_f = [];
arg.R2p_s = [];
arg.exchg = false; 

% substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% constant declarations
odims = size(ff);
nf = length(flip);
if arg.exchg
  L = 7;                        % free parameters: M0, ff, T1f, T1s, T2f, T2s, kfs
else
  L = 6;                        % free parameters: M0, ff, T1f, T1s, T2f, T2s
end

% if no mask specified, extrapolate to all voxels
if isempty(arg.mask)
  arg.mask = true(odims);
  N = prod(odims);
else
  N = numel(arg.mask(arg.mask));
end

% if no M0 specified, set to ones(odims) to yield normalized crb
if isempty(arg.M0)
  arg.M0 = ones(odims);
elseif norm(arg.M0 - ones(odims)) > N*eps
  warning('M0 other than unity not recommended: Fisher matrix may be scaled incorrectly.');
end

% if kfs specified, make sure exchg is on
% if no kfs specified, make sure exchg is off and set to zero
if isempty(arg.kfs)
  if arg.exchg
    warn('Exchange modeled but not specified? kfs/ksf set to zero.');
  end
  arg.kfs = zeros(odims);
elseif ~arg.exchg && nnz(arg.kfs) > 0
  warn('Since kfs set and nonzero, now modeling with exchange.');
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
M0 = masker(arg.M0, arg.mask);
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

% construct precision (inverse covariance) matrix
prec = Gdiag(div0(1, [sigp; sigm].^2));         % [D D]

% construct row gradient array
rgrad_p = NaN(N, L, nf);                        % [N L nf]
rgrad_m = NaN(N, L, nf);                        % [N L nf]
for i = 1:nf
  [rgrad_p(:,:,i), rgrad_m(:,:,i)] = dess_2comp_gradx(...
    M0, ff, T1f, T1s, T2f, T2s, kap, Dwf, Dws, flip(i), TR(i), TEp(i), TEm(i),...
    'kfs', kfs, 'ksf', ksf, 'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
end

% compute fisher matrix, voxel by voxel
F = NaN(N, L, L);                               % [N L L]
for n = 1:N
  tmp = cat(3, rgrad_p(n,:,:), rgrad_m(n,:,:)); % [1 L D]
  tmp = squeeze(tmp);                           % [L D]
  F(n,:,:) = real(tmp * prec * tmp');           % [L L]
end

% embed fisher matrices back into original mask
F = embed(F, arg.mask);                         % [(odims) L L]
end
