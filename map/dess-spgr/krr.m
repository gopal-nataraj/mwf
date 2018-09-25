  function [x, t] = krr(y, x0, nu, P, w, rff, train, dist, reg, dim, bool, mask)
%|function [x, t] = krr(y, x0, nu, P, w, rff, train, dist, reg, dim, bool, mask)
%|
%|  latent parameter estimation via kernel ridge regression
%|    approximates gaussian kernel via random fourier features
%|
%|  inputs
%|    y         [DV]            coil-combined image data
%|    x0        {1-L cell}      latent parameter initial estimates
%|                                any subset of {m0,ff,t1f,t1s,t2f,t2s}  
%|    nu        {N cell}        known parameters
%|    P         [1x1 struct]    scan parameters
%|                                1st field: data type (sp,de)
%|                                2nd field: scan parameter, as appropriate
%|     .*.tr    [S.*]               repetition times                                          ms
%|     .*.te    [S.* E.*]           echo times                                                ms
%|     .*.aex   [S.*]               nominal flip angle of excitation                          rad
%|    w         [D]             dataset weights
%|    rff       [1x1 struct]    random fourier features object
%|     .snr     [1]               estimate of max sig for unity m0
%|     .std     [1]               noise std dev in training data
%|     .len     [D+N]             kernel input length scales
%|     .c       [1]               global kernel length scale parameter 
%|     .H       [1]               embedding dimension       
%|     .K       [1]               number of training samples
%|    dist    	[1x1 struct]    sampling distribution object (ignored if x0 set!) 
%|     .x       {L cell}          latent parameter object
%|      .supp   [2]                 [lb ub] distribution support        
%|      .prior  {1}                 ('unif', 'logunif') distribution     
%|     .nu      {K cell}          known parameter object
%|      .supp   [2]                 [lb ub] distribution support    
%|    train     [1x1 struct]    krr training parameter object (if empty, will train)
%|     .mean.z  [H]               sample mean of feature maps
%|     .mean.x  [L]               sample mean of x
%|     .cov.zz  [H H]             sample auto-cov of feature maps
%|     .cov.xz  [L H]             sample cross-cov b/w x and feature maps
%|     .cov.xx  [L L]             sample auto-cov of latent parameters x
%|     .freq    [H D+N]           random 'frequency' vector
%|     .ph      [H]               random phase vector
%|    reg       [1]             krr regularization parameter
%|    dim       [1x1 struct]    object containing dimension info
%|    bool      [1x1 struct]    boolean variables
%|     .reset   false|true        reset random number generator during training
%|     .exchg   false|true        estimate exchange map                
%|     .mag.*   false|true        using magnitude (spgr,dess) data
%|     .chat    false|true        verbosity   
%|     .reset   false|true        reset random number generator during training
%|     .rfftst  false|true        show kernel approximation
%|     .nuclip  false|true        clip nu sampling distribution 
%|    mask      [(odims)]       binary object mask
%|  
%|  outputs
%|    x         {L cell}        latent object parameters
%|    t         [1]             (testing-only) run time
%|
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2016-09-01      original, based on krr_test.m
%|    1.2       2016-09-18      zero out nan values in m0 estimate
%|    2.1       2017-02-08      using different regularization for each parameter
%|    3.1       2018-01-22      rff.snr now controls m0 distribution sampling
%|                              now passing krr regularization parameter
%|                              added krr nu distribution clipping

% check that data is pure real
if ~isreal(y)
  error('Kernel regression with random fourier features for pure-real data only!');
end

% transfer prespecified x0 maps as known parameters
x = cell(dim.L,1);
known = false(dim.L,1);
for l=1:dim.L
  if ~isempty(x0{l})
    if bool.chat
      fprintf('Detected prespecified l=%u latent parameter; omitting from ML estimation...\n', l);
    end
    x{l} = masker(x0{l}, mask);
    nu{end+1} = x{l}; %#ok<AGROW>
    known(l) = true;
  else
    x{l} = NaN(dim.V,1);
  end
end

% check rff.len reflects possible length changes in nu
if length(rff.len) ~= (dim.D+length(nu))
  error('Length mismatch between rff.len and [y; nu]!?');
end

% testing inputs
tmp = reshape(y, [dim.D dim.V]);
tmp2 = transpose([nu{:}]);                                               % [N V]
q = [tmp; tmp2];                                                         % [D+N V]

% check if training needed
tmp = isfield(train, {'mean', 'cov', 'freq', 'ph'});
if ~all(tmp)
  if bool.chat
    fprintf('Sample statistics not detected: now training...');
  end
  
  % train to estimate sample statistics
  tic;
  train = krr_train(rff, dist, y, w, nu, P, dim, bool);
  tmp = toc;
  
  if bool.chat
    fprintf('done in %0.3f s.\n', tmp);
  end
end

% feature maps
tic;
z = rff_map(q, rff, w, train.freq, train.ph, dim);                      % [H V]

% kernel ridge regression
%   x{1} improperly normalized here
%   spurious values required to call rgradx_f(...)
tmp = bsxfun(@minus, z, train.mean.z); 
tmp = (train.cov.zz + reg*speye(rff.H)) \ tmp;                          % [H V]
for l = 1:dim.L
  if ~known(l)
    x{l} = col(train.cov.xz(l,:) * tmp);                                % [V]
    x{l} = bsxfun(@plus, x{l}, train.mean.x(l));                        % [V]
  end
end
t = toc;

% least-squares estimate of m0
% tmp = rgradx_f(x, nu, P, dim, bool);                                    % [DV LV sparse]
% Am = tmp(:,1:dim.L:end);                                                % [DV V sparse]
% W = spdiags(col(repmat(w, [1 dim.V])), 0, dim.D*dim.V, dim.D*dim.V);    % [DV DV sparse]
% num = Am'*W*y;                                                          % [V]
% den = (Am'*W*Am)*ones(dim.V,1);                                         % [V]
% x{1} = div_safe(num, den);                                              % [V]
% x{1}(isnan(x{1})) = 0;
end