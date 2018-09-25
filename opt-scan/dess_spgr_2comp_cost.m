  function cost = dess_spgr_2comp_cost(P, varargin)
%|function cost = dess_spgr_2comp_cost(P, varargin)
%|
%|  cost function for optimized scan design for 2-compartment parameter estimation
%|    options available for optimizing w.r.t. (relative) standard deviation, (r)std
%|    assumes fs and ksf are fixed via constraints
%|    assumes kap (flip scaling), Dwf/s (off-res mean), R2pf/s (off-res bw) all known
%|    if exchg on/off, then L=7/6 latent parameters and K=5/6 known parameters
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
%|                                if L=7, include exchange effects
%|                                if L=6, neglect exchange effects
%|    x         [1x1 struct]    latent object parameter distributions
%|                                1st arg: M0, ff, T1f, T1s, T2f, T2s, (kfs)          
%|                                2nd arg: one of the following:
%|     .*.prior [char]              def: /'gauss'/'unif'/'gauss'/'gauss'/'gauss'/'gauss'/'unif'/
%|     .*.mean  [1]                 def: /1/ /400/1000/20/80/0.00/              / / /ms/ms/ms/ms/kHz/
%|     .*.std   [1]                 def: 20% of corresponding means             / / /ms/ms/ms/ms/kHz/
%|     .*.minmax[2]                 def: / /[0.03 0.21]'/ / / / / /             / / /ms/ms/ms/ms/kHz/
%|     .*.nsamp [1]                 def: /1/5/3/3/3/3/1/
%|     .*.cfvar f|t                 optimize rel std (true) vs. std (false)     def: false unless wght>0
%|    nu        [1x1 struct]    known object parameter distributions
%|                                1st arg: kap, Dwf, Dws, R2pf, R2ps, (kfs)
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
%|    cost      [1]             estimate of expected cost
%|
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2016-04-05      original, bugs galore
%|    1.2       2016-04-08      bugs fixed, enough to match test script for singleton-grid case
%|    1.3       2016-04-15      minor changes in calling dess/spgr_2comp_crb(...)
%|    1.4       2016-08-16      changed format of P
%|    1.5       2017-05-18      now checks Fisher matrix conditioning before inversion

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

arg.x.kfs.prior = 'gauss';
arg.x.kfs.mean = 0;
arg.x.kfs.std = [];
arg.x.kfs.minmax = [];
arg.x.kfs.nsamp = 1;
arg.x.kfs.cfvar = [];

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

arg.L = length(arg.wght);
if arg.L==7
  % include exchange as unknown
  arg.exchg = true;
elseif arg.L==6
  % neglect exchange as known
  arg.exchg = false;
  arg.nu.kfs = arg.x.kfs;
  arg.x = rmfield(arg.x, 'kfs');
end

% to enumerate conveniently, switch from struct to cell
arg.x = struct2cell(arg.x);
arg.nu = struct2cell(arg.nu);
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

% set cost scaling based on tm_cmp
if arg.tm_cmp
  scale = sum([P.de.tr; P.sp.tr]);
else
  scale = 1;
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

% evaluate weighted trace of precision matrix, over large grid
wtprec = populate_array(@w_trace_prec,...
 arg.dims, arg);                % [(dims)]

% expected cost (coarsely) integrates wtprec over distribution
cost = scale * sum(col(wtprec .* joint));
end


%   function array = populate_array(fn, dims, arg)
% %|function array = populate_array(fn, dims, arg)
% %|
% %|  evaluates a function element-wise on an array
% %|    requires recursive calls to populate_array_helper(...)
% %|
% %|  input
% %|    fn        [func handle]   element-wise function operator
% %|                                assumes form fn(array, idx, arg) 
% %|                                evaluates on array at row vector position idx                              
% %|    dims      [1 ndims]       output array dimensions
% %|    arg       [struct]        obj containing add'l params needed in fn
% %|
% %|  output
% %|    array     [(dims)]        fn evaluated at every possible idx
% %|
% %|  version control
% %|    1.1       2016-04-07      original
% 
% % instantiate array
% array = NaN(dims);              % [(dims)]
% 
% % root call to recursive helper
% array = populate_array_helper(fn, array, dims, [], 1, arg);
% end
% 
% 
%   function array = populate_array_helper(fn, array, dims, idx, curr, arg)
% %|function array = populate_array_helper(fn, array, dims, idx, curr, arg)
% %|
% %|  helper function to populate_array(...)
% %|
% %|  input
% %|    fn        [func handle]   element-wise function operator
% %|                                assumes form fn(array, idx, arg) 
% %|                                evaluates on array at row vector position idx
% %|    array     [(dims)]        output array, progressively populated
% %|    dims      [1 0-ndims]     subset of output array dimensions
% %|    idx       [1 0-ndims]     subset of current index at which to eval fn
% %|    curr      [1]             which dim currently being iterated
% %|    arg       [struct]        obj containing add'l params needed in fn
% %|
% %|  version control
% %|    1.1       2016-04-07      original
% 
% % base case: evaluate fn, with idx of size [1 ndims]
% if isempty(dims)
%   array = fn(array, idx, arg);
% % recursive step: continue to build idx of size [1 <ndims]
% else
%   rest = dims(2:end);
%   for d = 1:dims(1)
%     if length(idx) < curr
%       idx = [idx d];
%     else
%       idx(curr) = d;
%     end
%     array = populate_array_helper(fn, array, rest, idx, curr+1, arg);
%   end
% end
% end

  function out = w_trace_prec(idx, arg)  
%|function out = w_trace_prec(idx, arg)
%| 
%|  evaluates weighted trace of precision matrix
%|
%|  input
%|    idx       [1 ndims]       index to select which arg.x.val and arg.nu.val
%|    arg       [struct]        obj containing params from parent fn
%|
%|  output
%|    out       [1]             weighted trace of precision matrix
%|
%|  version control
%|    1.1       2016-04-07      original
%|    1.2       2016-04-20      array handling removed
 
% separate indices as appropriate
idx_x  = idx(1, 1:arg.L);
idx_nu = idx(1, arg.L+1:end);

if arg.exchg
  % with exchg modeled, extract kfs from x{7}
  kfs_val = arg.x{7}.val(idx_x(7));
else
  % with exchg neglected, extract kfs from nu{6}
  kfs_val = arg.nu{6}.val(idx_nu(6));
end

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
  'kfs', kfs_val,...
  'R2p_f', arg.nu{4}.val(idx_nu(4)),...
  'R2p_s', arg.nu{5}.val(idx_nu(5)),...
  'exchg', arg.exchg);

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
  'kfs', kfs_val,...
  'R2p_f', arg.nu{4}.val(idx_nu(4)),...
  'R2p_s', arg.nu{5}.val(idx_nu(5)),...
  'exchg', arg.exchg);
    
% total information
F = squeeze(Fd + Fs);           % [L L]

% check condition number
if cond(F)>arg.cond
  error('Fisher matrix condition number in excess of %3e!?\nCheck initialization...',...
    arg.cond);
end

% assign weighted precision
out = trace(arg.W * (F \ (arg.W')));
end
