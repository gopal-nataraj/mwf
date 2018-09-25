  function grad = dess_spgr_2comp_costgrad(P, varargin)
%|function grad = dess_spgr_2comp_costgrad(P, varargin)
%|
%|  gradient of cost function for optimized scan design for 2-compartment parameter estimation
%|    options available for computing w.r.t. (relative) standard deviation, (r)std
%|    assumes fs fixed via constraint ff + fs = 1
%|    assumes kap (flip scaling), Dwf/s (off-res mean), R2pf/s (off-res bw) all known
%|    neglects exchange entirely
%|
%|  inputs
%|    P         [1x1 struct]    scan parameters
%|                                1st field: data type (sp,de)
%|                                2nd field: scan parameter, as appropriate
%|     .*.tr    [S.*]               repetition times                            ms
%|     .*.aex   [S.*]               nominal flip angle of excitation            rad
%|
%|  options
%|    wght      [L]             relative importance of estimating parameters    def: [0;1;0;0;0;0]
%|    x         [1x1 struct]    latent object parameter distributions
%|                                1st arg: M0, ff, T1f, T1s, T2f, T2s          
%|                                2nd arg: one of the following:
%|     .*.prior [char]              def: /'gauss'/'unif'/'gauss'/'gauss'/'gauss'/'gauss'/
%|     .*.mean  [1]                 def: /1/ /400/1000/20/80/                   / /ms/ms/ms/ms/kHz/
%|     .*.std   [1]                 def: 20% of corresponding means             / /ms/ms/ms/ms/kHz/
%|     .*.minmax[2]                 def: / /[0.03 0.21]'/ / / / /               / /ms/ms/ms/ms/kHz/
%|     .*.nsamp [1]                 def: /1/5/3/3/3/3/
%|     .*.cfvar f|t                 optimize rel std (true) vs. std (false)     def: false unless wght>0
%|    nu        [1x1 struct]    known object parameter distributions
%|                                1st arg: kap, Dwf, Dws, R2pf, R2ps
%|                                2nd arg: one of the following:
%|     .*.prior [char]              def: /'unif'/'gauss'/'gauss'/'gauss'/'gauss'/
%|     .*.mean  [1]                 def: / /0/0/0/0/                            / /kHz/kHz/kHz/kHz/
%|     .*.std   [1]                 def: /20%/eps/eps/eps/eps/                  / /kHz/kHz/kHz/kHz/
%|     .*.minmax[2]                 def: /[0.9 1.1]'/ / / / /                   / /kHz/kHz/kHz/kHz/
%|     .*.nsamp [1]                 def: /3/1/1/1/1/
%|     .*.cfvar f|t                 optimize rel std (true) vs. std (false)     def: false
%|    TE        [1x1 struct]    echo time objects (typically all same)
%|     .p       [S.de]             dess 'defocusing' echo time                  def: 4.67ms
%|     .m       [S.de]             dess 'refocusing' echo time                  def: 4.67ms
%|     .s       [S.sp]             spgr 'defocusing' echo time                  def: 4.67ms
%|    sig       [1x1 struct]    unit-M0 noise std deviations
%|     .p       [S.de]             dess 'defocusing' echo noise std dev.        def: 3.8607e-4
%|     .m       [S.de]             dess 'refocusing' echo noise std dev.        def: 3.8607e-4
%|     .s       [S.sp]             spgr 'defocusing' echo noise std dev.        def: 3.8607e-4
%|    cond      [1]             max permissible Fisher matrix conditionality    def: 10^20
%|    tm_cmp    f|t             overall cost time compensation off|on           def: false
%|
%|  outputs
%|    grad      [1x1 struct]    gradients of cost at P
%|                                1st field: data type (sp,de)
%|                                2nd field: scan parameter, as appropriate
%|     .*.tr    [S.*]               repetition times                            ms
%|     .*.aex   [S.*]               nominal flip angle of excitation            rad
%|
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2016-04-20      original
%|    1.2       2016-04-26      corrected to handle joint distributions larger than size [1]
%|    1.3       2016-08-16      changed format of P and grad
%|    1.4       2017-05-18      now checks Fisher matrix conditioning before inversion

% constant declarations
S.de = length(P.de.aex);
S.sp = length(P.sp.aex);

% error checks on TRd, TRs
if length(P.de.tr) ~= S.de
  error('flipd and TRd of unequal length!');
elseif length(P.sp.tr) ~= S.sp
  error('flips and TRs of unequal length!');
end

% default x values
arg.x.M0.prior = 'gauss';
arg.x.M0.mean = 1;
arg.x.M0.std = [];
arg.x.M0.minmax = [];
arg.x.M0.nsamp = 1;
arg.x.M0.cfvar = [];

arg.x.ff.prior = 'unif';
arg.x.ff.mean = [];
arg.x.ff.std = [];
arg.x.ff.minmax = [0.03 0.21]';
arg.x.ff.nsamp = 5;
arg.x.ff.cfvar = [];

arg.x.T1f.prior = 'gauss';
arg.x.T1f.mean = 400;
arg.x.T1f.std = [];
arg.x.T1f.minmax = [];
arg.x.T1f.nsamp = 3;
arg.x.T1f.cfvar = [];

arg.x.T1s.prior = 'gauss';
arg.x.T1s.mean = 1000;
arg.x.T1s.std = [];
arg.x.T1s.minmax = [];
arg.x.T1s.nsamp = 3;
arg.x.T1s.cfvar = [];

arg.x.T2f.prior = 'gauss';
arg.x.T2f.mean = 20;
arg.x.T2f.std = [];
arg.x.T2f.minmax = [];
arg.x.T2f.nsamp = 3;
arg.x.T2f.cfvar = [];

arg.x.T2s.prior = 'gauss';
arg.x.T2s.mean = 80;
arg.x.T2s.std = [];
arg.x.T2s.minmax = [];
arg.x.T2s.nsamp = 3;
arg.x.T2s.cfvar = [];

% default nu values
arg.nu.kap.prior = 'unif';
arg.nu.kap.mean = [];
arg.nu.kap.std = [];
arg.nu.kap.minmax = [0.9 1.1]';
arg.nu.kap.nsamp = 3;
arg.nu.kap.cfvar = false;

arg.nu.Dwf.prior = 'gauss';
arg.nu.Dwf.mean = 0;
arg.nu.Dwf.std = [];
arg.nu.Dwf.minmax = [];
arg.nu.Dwf.nsamp = 1;
arg.nu.Dwf.cfvar = false;

arg.nu.Dws.prior = 'gauss';
arg.nu.Dws.mean = 0;
arg.nu.Dws.std = [];
arg.nu.Dws.minmax = [];
arg.nu.Dws.nsamp = 1;
arg.nu.Dws.cfvar = false;

arg.nu.R2pf.prior = 'gauss';
arg.nu.R2pf.mean = 0;
arg.nu.R2pf.std = [];
arg.nu.R2pf.minmax = [];
arg.nu.R2pf.nsamp = 1;
arg.nu.R2pf.cfvar = false;

arg.nu.R2ps.prior = 'gauss';
arg.nu.R2ps.mean = 0;
arg.nu.R2ps.std = [];
arg.nu.R2ps.minmax = [];
arg.nu.R2ps.nsamp = 1;
arg.nu.R2ps.cfvar = false;

% other default values
arg.wght = [0,1,0,0,0,0]';
arg.TE.p = 4.67 * ones(S.de,1);  % ms
arg.TE.m = 4.67 * ones(S.de,1);  % ms
arg.TE.s = 4.67 * ones(S.sp,1);  % ms
arg.sig.p = 3.8607e-4 * ones(S.de,1);
arg.sig.m = 3.8607e-4 * ones(S.de,1);
arg.sig.s = 3.8607e-4 * ones(S.sp,1);
arg.cond = 10^20;
arg.tm_cmp = false;

% substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% error checks on weights
if sum(arg.wght<0) > 0
  error('negative weights!');
elseif norm(arg.wght,1) ~= 1
  tmp = norm(arg.wght,1);
  arg.wght = div0(arg.wght, tmp);
  warning('wghts renormalized by factor of %0.2f', tmp);
end

% to enumerate conveniently, switch from struct to cell
arg.x = struct2cell(arg.x);
arg.nu = struct2cell(arg.nu);
arg.L = length(arg.x);
arg.K = length(arg.nu);

% error checks and default values for x
for l = 1:arg.L
  % use coeff of variation only if corresponding wght>0 set
  if isempty(arg.x{l}.cfvar)
    if arg.wght(l) > 0
      arg.x{l}.cfvar = true;
    else
      arg.x{l}.cfvar = false;
    end
  end
  % set mean and std dev of uniform prior, using minmax
  if strcmp(arg.x{l}.prior, 'unif')
    if isempty(arg.x{l}.minmax)
      error('for l=%u latent param: uniform dist requires minmax set.');
    elseif isempty(arg.x{l}.mean)
      arg.x{l}.mean = mean(arg.x{l}.minmax);
    end
    if isempty(arg.x{l}.std)
      arg.x{l}.std = 0.2*arg.x{l}.mean;
    end
    % reset minmax to mean if only one sample
    if arg.x{l}.nsamp == 1
      arg.x{l}.minmax = arg.x{l}.mean * [1 1]';
    end
  % set minmax and std dev of gaussian prior, using mean
  elseif strcmp(arg.x{l}.prior, 'gauss')
    if isempty(arg.x{l}.mean)
      error('for l=%u latent param: gaussian dist requires mean set.');
    elseif isempty(arg.x{l}.std)
      arg.x{l}.std = max(abs(0.2*arg.x{l}.mean), eps);
    end
    if isempty(arg.x{l}.minmax)
      % sample over 1-3 std dev, depending on nsamp
      nsd = min(floor(div0(arg.x{l}.nsamp,2)), 3);
      arg.x{l}.minmax = arg.x{l}.mean + arg.x{l}.std * [-nsd nsd]';
    end
  else
    error('For l=%u latent param: invalid distribution!', l);
  end
end

% error checks and default values for nu
for k = 1:arg.K
  % set mean and std dev of uniform prior, using minmax
  if strcmp(arg.nu{k}.prior, 'unif')
    if isempty(arg.nu{k}.minmax)
      error('for k=%u known param: uniform dist requires minmax set.');
    elseif isempty(arg.nu{k}.mean)
      arg.nu{k}.mean = mean(arg.nu{k}.minmax);
    end
    if isempty(arg.nu{k}.std)
      arg.nu{k}.std = 0.2*arg.nu{k}.mean;
    end
    % reset minmax to mean if only one sample
    if arg.nu{k}.nsamp == 1
      arg.nu{k}.minmax = arg.nu{k}.mean * [1 1]';
    end
  % set minmax and std dev of gaussian prior, using mean
  elseif strcmp(arg.nu{k}.prior, 'gauss')
    if isempty(arg.nu{k}.mean)
      error('for k=%u known param: gaussian dist requires mean set.');
    elseif isempty(arg.nu{k}.std)
      arg.nu{k}.std = max(abs(0.2*arg.nu{k}.mean), eps);
    end
    if isempty(arg.nu{k}.minmax)
      % sample over 1-3 std dev, depending on nsamp
      nsd = min(floor(div0(arg.nu{k}.nsamp,2)), 3);
      arg.nu{k}.minmax = arg.nu{k}.mean + arg.nu{k}.std * [-nsd nsd]';
    end
  else
    error('For k=%u known param: invalid distribution!', k);
  end
end

% error checks on std dev lengths
if length(arg.sig.p) ~= S.de || length(arg.sig.m) ~= S.de
  error('sigp and sigm must be of same length as flipd!');
elseif length(arg.sig.s) ~= S.sp
  error('sigs must be of same length as flips!');
end

% if using time compensation, set scaling
if arg.tm_cmp
  scale = sum([P.de.tr; P.sp.tr]);
end

% design x prior
for l = 1:arg.L
  % make evaluation points
  arg.x{l}.val = col(linspace(...
    arg.x{l}.minmax(1), arg.x{l}.minmax(2), arg.x{l}.nsamp));
        
  % make prior distributions
  if strcmp(arg.x{l}.prior, 'unif')
    arg.x{l}.dist = div0(1,arg.x{l}.nsamp) * ones(arg.x{l}.nsamp,1);
  elseif strcmp(arg.x{l}.prior, 'gauss')
    arg.x{l}.dist = normpdf(arg.x{l}.val, arg.x{l}.mean, arg.x{l}.std);
  end

  % compensate, if using coefficients of variation
  % trick: weight by at most 1/eps
  if arg.x{l}.cfvar
    arg.x{l}.dist = div0(arg.x{l}.dist, max(abs(arg.x{l}.val).^2, eps));
  end 

  % normalize distributions
  % trick, cont'd: essentially a delta function if val contains zero
  arg.x{l}.dist = div0(arg.x{l}.dist, norm(arg.x{l}.dist,1));
end

% design nu prior
for k = 1:arg.K
  % make evaluation points
  arg.nu{k}.val = col(linspace(...
    arg.nu{k}.minmax(1), arg.nu{k}.minmax(2), arg.nu{k}.nsamp));
    
  % make prior distributions
  if strcmp(arg.nu{k}.prior, 'unif')
    arg.nu{k}.dist = div0(1,arg.nu{k}.nsamp) * ones(arg.nu{k}.nsamp,1);
  elseif strcmp(arg.nu{k}.prior, 'gauss')
    arg.nu{k}.dist = normpdf(arg.nu{k}.val, arg.nu{k}.mean, arg.nu{k}.std);
  end
    
  % compensate, if using coefficients of variation
  % trick: weight by at most 1/eps
  if arg.nu{k}.cfvar
    arg.nu{k}.dist = div0(arg.nu{k}.dist, max(abs(arg.nu{k}.val).^2, eps));
  end
    
  % normalize distribution
  % trick, cont'd: essentially a delta function if val contains zero
  arg.nu{k}.dist = div0(arg.nu{k}.dist, norm(arg.nu{k}.dist,1));
end

% construct joint distribution from marginals
% note: dims = [x{1}.nsamp...x{L}.nsamp} nu{1}.nsamp...n{K}.nsamp]
vecs = cell(arg.L+arg.K, 1);
arg.dims = NaN(1, arg.L+arg.K);
for l = 1:arg.L
  vecs{l} = arg.x{l}.dist;
  arg.dims(l) = arg.x{l}.nsamp;
end
for k = 1:arg.K
  vecs{arg.L+k} = arg.nu{k}.dist;
  arg.dims(arg.L+k) = arg.nu{k}.nsamp;
end
joint = outer_prod(vecs);       % [(dims)]

% normalization check on joint distribution
if norm(sum(col(joint)) - 1) > eps*numel(joint)
  error('Joint distribution not normalized?');
end

% construct diagonal weighting matrix
% since small, use full matrix to later use backslash
arg.W = diag(arg.wght);         % [L L]

% save input P as a field of arg for porting
arg.P = P;

% evaluate gradient struct array, over large grid
% each element is a grad struct
elem.flipd = NaN(S.de,1);
elem.flips = NaN(S.sp,1);
elem.TRd = NaN(S.de,1);
elem.TRs = NaN(S.sp,1);
gradArr = pop_struct_array(@cost_grad,...
  elem, arg.dims, arg);         % [(dims) struct]

% extract individual parameters from gradient struct
flipdArray  = [gradArr.flipd];  % [S.de prod(odims)]
flipsArray  = [gradArr.flips];  % [S.sp prod(odims)]
TRdArray    = [gradArr.TRd];    % [S.de prod(odims)]
TRsArray    = [gradArr.TRs];    % [S.sp prod(odims)]

% compute expected gradient
grad.de.aex = ...               % [S.de]
  sum(flipdArray  .* repmat(joint(:).', [S.de 1]), 2);
grad.sp.aex = ...               % [S.sp]
  sum(flipsArray  .* repmat(joint(:).', [S.sp 1]), 2);
grad.de.tr = ...                % [S.de]
  sum(TRdArray    .* repmat(joint(:).', [S.de 1]), 2);
grad.sp.tr = ...                % [S.sp]
  sum(TRsArray    .* repmat(joint(:).', [S.sp 1]), 2);

% if scaling on, modify gradient via product rule
if arg.tm_cmp
  % evaluate cost w/ same options, except no time compensation
  costArg = cell_rep(varargin, 'tm_cmp', false);
  cost = dess_spgr_2comp_cost(arg.P, costArg);

  % product rule
  grad.de.aex = scale * grad.de.aex;
  grad.sp.aex = scale * grad.sp.aex;
  grad.de.tr  = scale * grad.de.tr + cost*ones(S.de,1);
  grad.sp.tr  = scale * grad.sp.tr + cost*ones(S.sp,1);
end
end
  
  
  function gradP = cost_grad(idx, arg)
%|function gradP = cost_grad(idx, arg)
%|
%|  evaluates cost function gradient
%|
%|  input
%|    idx       [1 ndims]       index to select which arg.x.val and arg.nu.val
%|    arg       [struct]        obj containing params from parent fn
%|
%|  output
%|    gradP     [1x1 struct]    gradients of cost at P
%|      .flipd  [S.de]            gradients w.r.t. flipd
%|      .flips  [S.sp]            gradients w.r.t. flips
%|      .TRd    [S.de]            gradients w.r.t. TRd
%|      .TRs    [S.sp]            gradients w.r.t. TRs
%|
%|  version control
%|    1.1       2016-04-20      original
%|    1.2       2016-08-16      changed format of P

% instantiate gradient object
gradP.flipd   = NaN(length(arg.P.de.aex), 1);
gradP.flips   = NaN(length(arg.P.sp.aex), 1);
gradP.TRd     = NaN(length(arg.P.de.aex), 1);
gradP.TRs     = NaN(length(arg.P.sp.aex), 1);

% separate indices as appropriate
idx_x  = idx(1, 1:arg.L);
idx_nu = idx(1, arg.L+1:end);

% dess 2-component fisher information
Fd = dess_2comp_crb(...         % [1 L L]
  arg.x{2}.val(idx_x(2)),...
  arg.x{3}.val(idx_x(3)),...
  arg.x{4}.val(idx_x(4)),...
  arg.x{5}.val(idx_x(5)),...
  arg.x{6}.val(idx_x(6)),...
  arg.nu{1}.val(idx_nu(1)),...
  arg.nu{2}.val(idx_nu(2)),...
  arg.nu{3}.val(idx_nu(3)),...
  arg.P.de.aex,...
  arg.P.de.tr,...
  arg.TE.p,...
  arg.TE.m,...
  arg.sig.p,...
  arg.sig.m,...
  'M0', arg.x{1}.val(idx_x(1)),...
  'R2p_f', arg.nu{4}.val(idx_nu(4)),...
  'R2p_s', arg.nu{5}.val(idx_nu(5)),...
  'exchg', false);

% spgr 2-component fisher information
Fs = spgr_2comp_crb(...         % [1 L L]
  arg.x{2}.val(idx_x(2)),...
  arg.x{3}.val(idx_x(3)),...
  arg.x{4}.val(idx_x(4)),...
  arg.x{5}.val(idx_x(5)),...
  arg.x{6}.val(idx_x(6)),...
  arg.nu{1}.val(idx_nu(1)),...
  arg.nu{2}.val(idx_nu(2)),...
  arg.nu{3}.val(idx_nu(3)),...
  arg.P.sp.aex,...
  arg.P.sp.tr,...
  arg.TE.s,...
  arg.sig.s,...
  'M0', arg.x{1}.val(idx_x(1)),...
  'R2p_f', arg.nu{4}.val(idx_nu(4)),...
  'R2p_s', arg.nu{5}.val(idx_nu(5)),...
  'exchg', false);

% total information
F = squeeze(Fd + Fs);           % [L L]

% check condition number
if cond(F)>arg.cond
  error('Fisher matrix condition number in excess of %3e!?\nCheck initialization...',...
    arg.cond);
end

% intermediate variable
M = arg.W / F;

% compute gradient dataset-by-dataset
for ifd = 1:length(arg.P.de.aex)
  % dess signal gradient w.r.t. x
  [sp_gradx, sm_gradx] = dess_2comp_gradx(...
    arg.x{1}.val(idx_x(1)),...
    arg.x{2}.val(idx_x(2)),...
    arg.x{3}.val(idx_x(3)),...
    arg.x{4}.val(idx_x(4)),...
    arg.x{5}.val(idx_x(5)),...
    arg.x{6}.val(idx_x(6)),...
    arg.nu{1}.val(idx_nu(1)),...
    arg.nu{2}.val(idx_nu(2)),...
    arg.nu{3}.val(idx_nu(3)),...
    arg.P.de.aex(ifd),...
    arg.P.de.tr(ifd),...
    arg.TE.p(ifd),...
    arg.TE.m(ifd),...
    'R2p_f', arg.nu{4}.val(idx_nu(4)),...
    'R2p_s', arg.nu{5}.val(idx_nu(5)),...
    'exchg', false);            % [1 L]
  
  % dess signal mixed gradients w.r.t. x, p
  [sp_gradx_gradp, sm_gradx_gradp] = dess_2comp_gradx_gradp(...
    arg.x{1}.val(idx_x(1)),...
    arg.x{2}.val(idx_x(2)),...
    arg.x{3}.val(idx_x(3)),...
    arg.x{4}.val(idx_x(4)),...
    arg.x{5}.val(idx_x(5)),...
    arg.x{6}.val(idx_x(6)),...
    arg.nu{1}.val(idx_nu(1)),...
    arg.nu{2}.val(idx_nu(2)),...
    arg.nu{3}.val(idx_nu(3)),...
    arg.P.de.aex(ifd),...
    arg.P.de.tr(ifd),...
    arg.TE.p(ifd),...
    arg.TE.m(ifd),...
    'R2p_f', arg.nu{4}.val(idx_nu(4)),...
    'R2p_s', arg.nu{5}.val(idx_nu(5)),...
    'exchg', false);            % [1 L P]
  
  % intermediate variables
  vp = col(sp_gradx);           % [L 1]
  vpdot_fl = ...  
    col(sp_gradx_gradp(:,:,1)); % [L 1]
  vpdot_tr = ...
    col(sp_gradx_gradp(:,:,2)); % [L 1]
  
  vm = col(sm_gradx);           % [L 1]
  vmdot_fl = ...
    col(sm_gradx_gradp(:,:,1)); % [L 1]
  vmdot_tr = ...
    col(sm_gradx_gradp(:,:,2)); % [L 1]
  
  % cost gradient assigments
  tmp.p = div0(vp*vpdot_fl', arg.sig.p(ifd)^2);
  tmp.m = div0(vm*vmdot_fl', arg.sig.m(ifd)^2);
  gradP.flipd(ifd) = ...        % [1]
    -2 * real(trace(M * (tmp.p + tmp.m) * M'));
  
  tmp.p = div0(vp*vpdot_tr', arg.sig.p(ifd)^2);
  tmp.m = div0(vm*vmdot_tr', arg.sig.m(ifd)^2);
  gradP.TRd(ifd) = ...          % [1]
    -2 * real(trace(M * (tmp.p + tmp.m) * M'));
end
for ifs = 1:length(arg.P.sp.aex)
  % spgr signal gradient w.r.t. x
  [ss_gradx] = spgr_2comp_gradx(...
    arg.x{1}.val(idx_x(1)),...
    arg.x{2}.val(idx_x(2)),...
    arg.x{3}.val(idx_x(3)),...
    arg.x{4}.val(idx_x(4)),...
    arg.x{5}.val(idx_x(5)),...
    arg.x{6}.val(idx_x(6)),...
    arg.nu{1}.val(idx_nu(1)),...
    arg.nu{2}.val(idx_nu(2)),...
    arg.nu{3}.val(idx_nu(3)),...
    arg.P.sp.aex(ifs),...
    arg.P.sp.tr(ifs),...
    arg.TE.s(ifs),...
    'R2p_f', arg.nu{4}.val(idx_nu(4)),...
    'R2p_s', arg.nu{5}.val(idx_nu(5)),...
    'exchg', false);            % [1 L]
  
  % spgr signal mixed gradients w.r.t. x, p
  [ss_gradx_gradp] = spgr_2comp_gradx_gradp(...
    arg.x{1}.val(idx_x(1)),...
    arg.x{2}.val(idx_x(2)),...
    arg.x{3}.val(idx_x(3)),...
    arg.x{4}.val(idx_x(4)),...
    arg.x{5}.val(idx_x(5)),...
    arg.x{6}.val(idx_x(6)),...
    arg.nu{1}.val(idx_nu(1)),...
    arg.nu{2}.val(idx_nu(2)),...
    arg.nu{3}.val(idx_nu(3)),...
    arg.P.sp.aex(ifs),...
    arg.P.sp.tr(ifs),...
    arg.TE.s(ifs),...
    'R2p_f', arg.nu{4}.val(idx_nu(4)),...
    'R2p_s', arg.nu{5}.val(idx_nu(5)),...
    'exchg', false);            % [1 L P]
    
  % intermediate variables
  vs = col(ss_gradx);           % [L 1]
  vsdot_fl = ...
    col(ss_gradx_gradp(:,:,1)); % [L 1]
  vsdot_tr = ...
    col(ss_gradx_gradp(:,:,2)); % [L 1]
  
  % cost gradient assignments
  tmp.s = div0(vs*vsdot_fl', arg.sig.s(ifs)^2);
  gradP.flips(ifs) = ...        % [1]
    -2 * real(trace(M * tmp.s * M'));
  
  tmp.s = div0(vs*vsdot_tr', arg.sig.s(ifs)^2);
  gradP.TRs(ifs) = ...          % [1]
    -2 * real(trace(M * tmp.s * M'));
end
end
