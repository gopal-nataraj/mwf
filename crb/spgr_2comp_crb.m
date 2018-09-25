  function F = spgr_2comp_crb(ff, T1f, T1s, T2f, T2s,...
      kap, Dwf, Dws, flip, TR, TE, sig, varargin)
%|function F = spgr_2comp_crb(ff, T1f, T1s, T2f, T2s,...
%|    kap, Dwf, Dws, flip, TR, TE, sig, varargin)
%|
%|  2-compartment spgr signal fisher information matrix
%|    assumes fs and ksf are fixed via constraints
%|    assumes Dwf, Dws, kap (flip scaling) all known
%|    if exchg on/off, then L=7/6 total parameters
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
%|    TE        [nf]            echo time                                       ms
%|    sig       [nf]            unit-M0 noise standard deviation          
%|
%|  options
%|    mask      [(odims)]       object mask                                     def: true(odims)
%|    M0        [(odims)]       spin density (recommend not changing!)          def: ones(odims)
%|    fs        [(odims)]       slow component spin fraction                    def: 1-ff
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
%|    1.1       2016-03-31      original
%|    1.2       2016-04-15      modified to work with spgr_2comp_gradx(...)
%|    1.3       2017-05-25      now works properly in special case of nf==1
%|    1.4       2018-01-21      now omits fs as an independent variable (forcing fs := 1-ff)

% dafault values
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
% kfs = masker(arg.kfs, arg.mask);
% ksf = masker(arg.ksf, arg.mask);
kap = masker(kap, arg.mask);
Dwf = masker(Dwf, arg.mask);
Dws = masker(Dws, arg.mask);
R2pf = masker(arg.R2p_f, arg.mask);
R2ps = masker(arg.R2p_s, arg.mask);

% construct precision (inverse covariance) matrix
prec = Gdiag(div0(1, sig.^2));                  % [D D]

% construct row gradient array with analytical derivative
rgrad_s = NaN(N, L, nf);                        % [N L D]
for i = 1:nf
  if arg.exchg
    error('todo: gradient of signal with physical exchange');
  else
    rgrad_s(:,:,i) = spgr_2comp_gradx(...
      M0, ff, T1f, T1s, T2f, T2s, kap, Dwf, Dws, flip(i), TR(i), TE(i),...
      'R2p_f', R2pf, 'R2p_s', R2ps, 'exchg', arg.exchg);
  end
end

% compute fisher matrix, voxel-by-voxel
F = NaN(N, L, L);
for n = 1:N
  tmp = squeeze(rgrad_s(n,:,:));                % [L D], or [1 L] if D==1
  if nf==1
    tmp = col(tmp);                             % [L 1] 
  end
  F(n,:,:) = real(tmp * prec * tmp');           % [L L]
end

% embed fisher matrices back into original mask
F = embed(F, arg.mask);                         % [(odims) L L]
end
