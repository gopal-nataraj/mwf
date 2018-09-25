  function [x, t] = mri_multicomp_map(y, P, varargin)
%|function [x, t] = mri_multicomp_map(y, P, varargin)
%|
%|  mri parameter estimation from multi-compartmental models
%|    handles coil-combined image data from the following pulse sequences:
%|      spoiled gradient-recalled echo
%|        assumes perfect spoiling
%|      dual-echo steady state
%|        neglects diffusion effects
%|    handles only two components presently
%|    permits b1, b0, and r2p variation
%|    neglects physical exchange
%|
%|  inputs
%|    y         [1x1 struct]    coil-combined image data
%|     .sp      [(odims) S.sp 1]  spoiled gradient-recalled echo
%|     .de      [(odims) S.de 2]  dual-echo steady state 
%|    P         [1x1 struct]    scan parameters
%|                                1st field: data type (sp,de)
%|                                2nd field: scan parameter, as appropriate
%|     .*.tr    [S.*]               repetition times                                          ms
%|     .*.te    [S.* E.*]           echo times                                                ms
%|     .*.aex   [S.*]               nominal flip angle of excitation                          rad
%|
%|  options
%|    mask      [1x1 struct]    binary object masks
%|     .disp    [(odims)]         over which to display images          def: true(odims)
%|     .est     [(odims)]         over which to perform estimation      def: imdilate(mask.disp)
%|     .noise   [(odims)]         over which to estimate noise std dev  def: ~mask.est
%|    nu        [1x1 struct]    known parameters
%|     .kap     [(odims)]         flip-angle scale map                  def: ones(odims)            
%|     .b0f     [(odims)]         fast off-resonance map                def: zeros(odims)     kHz
%|     .b0s     [(odims)]         slow off-resonance map                def: zeros(odims)     kHz
%|     .r2pf    [(odims)]         fast broadening linewidth map         def: zeros(odims)     kHz
%|     .r2ps    [(odims)]         slow broadening linewidth map         def: zeros(odims)     kHz
%|    wght      [1x1 struct]    dataset weights
%|     .sp      [S.sp 1]          spoiled GRE weights                   def: ones(S.sp,1)
%|     .de      [S.de 2]          dual echo steady-state weights        def: ones(S.de,1)
%|    thresh    [1]             frac of max sig deciding 'background' 	def: 0.05
%|    x0        [1x1 struct]    latent parameter initial estimates      def: from init
%|     .m0      [(odims)]         spin density
%|     .ff      [(odims)]         fast fraction
%|     .t1f     [(odims)]         fast spin-lattice relaxation time                           ms
%|     .t1s     [(odims)]         slow spin-lattice relaxation time                           ms
%|     .t2f     [(odims)]         fast spin-spin relaxation time                              ms
%|     .t2s     [(odims)]         slow spin-spin relaxation time                              ms
%|    meth      [1x1 struct]    estimation methods                  
%|     .init    {1}               initialization                        def: 'krr'
%|     .iter    {1}               iterative local optimization          def: 'pgpm' 
%|    dist.x    [1x1 struct]    latent parameter sampling distribution object (ignored if x0 set!)  
%|                                2nd field: latent parameter (m0,ff,t1f,t1s,t2f,t2s,(kfs))
%|     .*.supp  [2]               [lb ub] distribution support          def: see below
%|     .*.nsamp [1]               (vpm) number of dict samples          def: see below
%|     .*.prior {1}               ('unif', 'logunif') distribution      def: see below 
%|    dist.nu   [1x1 struct]    known parameter sampling distribution object
%|                                2nd field: known parameter (kap,b0f,b0s,r2pf,r2ps)
%|     .*.supp  [2]               [lb ub] distribution support          def: see below    
%|    kmean     [1x1 struct]    (vpm,pgpm) kmeans object to pool x0/nu
%|     .C       [1]               number of kmeans clusters             def: 10
%|     .opt     {1 2*nopt}        optional arguments to kmeans          def: see below
%|    rff       [1x1 struct]    (krr) random fourier features object
%|     .snr     [1]               estimate of max sig for unity m0      def: 0.1
%|     .std     [1]               noise std dev in training data        def: est from noise 
%|     .len     [D+N]             kernel input length scales            def: from data
%|     .c       [1]               global kernel length scale parameter  def: 2^0
%|     .H       [1]               embedding dimension                   def: 10^4
%|     .K       [1]               number of training samples            def: 10^6
%|    train     [1x1 struct]    (krr) training parameter object         def: from training
%|     .mean.z  [H]               sample mean of feature maps
%|     .mean.x  [L]               sample mean of x
%|     .cov.zz  [H H]             sample auto-cov of feature maps
%|     .cov.xz  [L H]             sample cross-cov b/w x and feature maps
%|     .cov.xx  [L L]             sample auto-cov of latent parameters x
%|     .freq    [H D+N]           random 'frequency' vector
%|     .ph      [H]               random phase vector
%|    inv       [1x1 struct]    matrix inversion reg strength object
%|     .krr     [1]               kernel ridge regression               def: 10^-8
%|     .lm      [1]               levenberg-marquardt                   def: 10^-8
%|    precon    [1x1 struct]    preconditioning parameter object        
%|     .update  [1]               number of iter to update precon       def: 5
%|     .hessreg [1]               hessian regularization parameter      def: 10^-10
%|    line      [1x1 struct]    backtracking line search object
%|     .step0   [1]               initial step size                     def: 1
%|     .prop    [1]               step size reduction factor            def: 0.5
%|    boxcon    [1x1 struct]    [lb ub] box constraints for iter opt
%|                                1st field: latent parameter (m0,ff,t1f,t1s,t2f,t2s,(kfs))
%|    reg       [1x1 struct]    edge-preserving (hyper3) regularizers (see Reg1.m)   
%|                                1st field: latent parameter 
%|     .*.pot   {1 npotarg}       potential function arguments          def: {'hyper3', delta}      
%|     .*.beta  [1]               strength                              def: scales with D
%|    stop      [1x1 struct]    iteration stopping criteria
%|     .iter    [1]               maximum number of iterations          def: 100
%|     .wghtx   [L]               provides weighting in weighted norm   def: from boxcon
%|     .tolx    [1]               compare vs wnorm(x-xprev)/wnorm(x)    def: 10^-7
%|    bool      [1x1 struct]    boolean variables
%|     .exchg   false|true        estimate exchange map                 def: false
%|     .mag.*   false|true        using magnitude (sp,de) data          def: true
%|     .chat    false|true        verbosity                             def: false
%|     .norm    false|true        normalize data for scale-invar reg    def: true
%|     .reset   false|true        (krr) reset rng while sampling        def: true
%|     .rfftst  false|true        (krr) show kernel approximation       def: false
%|     .nuclip  false|true        (krr) clip nu sampling distribtion    def: true
%|     .reg     false|true        use regularization                    def: false
%|     .precon  false|true        use preconditioner (hessian approx)   def: true
%|     .disp    false|true        show image updates                    def: false
%|    disp      [1x1 struct]    display ranges for image updates       
%|     .m0      [2]               abs(m0) display range                 def: [0 2*median(abs(m0))]
%|     .ff      [2]               ff display range                      def: [0 0.4]
%|     .t1f     [2]               t1f display range                     def: [0 1000]         ms
%|     .t1s     [2]               t1s display range                     def: [0 2000]         ms
%|     .t2f     [2]               t2f display range                     def: [0 50]           ms
%|     .t2s     [2]               t2s display range                     def: [0 200]          ms
%|
%|  outputs
%|    x         [1x1 struct]    object estimates
%|     .init.*  [(odims)]         initial estimates
%|     .iter.*  [(odims)]         iterative estimates, if computed
%|                                2nd field: latent parameter (m0,t1,t2(,inveff))
%|    t.*       [1x1 struct]    run times (init,iter)
%|
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2016-08-12      skeleton code
%|    1.2       2016-08-18      ml estimation working in simulation
%|    2.1       2016-08-21      added simulated annealing init option
%|    3.1       2016-09-01      added kernel regresssion init option
%|    3.2       2016-09-13      added automatic noise std dev estimation
%|    3.3       2016-09-16      make regularization optional
%|    3.4       2016-10-07      now handles case when x0 given fully
%|    3.5       2016-10-10      krr now samples m0 dist based on data
%|    4.1       2016-10-12      added precon gradient proj meth iterative option
%|    4.2       2016-10-20      changed default range of ff to [-0.1 0.4]
%|    5.1       2016-10-26      now using signal model hessians
%|    6.1       2017-02-08      using different regularization for each parameter
%|    6.2       2017-05-18      now sets dim.odims correctly in case of 1 spgr scan
%|    7.1       2018-01-22      rff.snr now controls m0 distribution sampling
%|                              added krr regularization strengh option
%|                              switched output format from e.g. x.m0.init to x.init.m0
%|                              added krr nu distribution clipping
%|                              added method-of-moments init option for spgr/dess
%|                              pgpm now loops over tissue clusters for separate line searches

% object dimensions
tmp = isfield(y, {'sp','de'});
if sum(tmp)==0
  error('Detected no data?!');
end
for i = 1:length(tmp)
  if tmp(i)
    switch i
      case 1
        tmp2 = size(y.sp);
        dim.odims = tmp2(1:end-1);
      case 2
        tmp2 = size(y.de);
        dim.odims = tmp2(1:end-2);
    end
  end
end

% create empty arrays for missing fields
for i = 1:length(tmp)
  if ~tmp(i)
    switch i
      case 1, y.sp = zeros([dim.odims 0 1]);
      case 2, y.de = zeros([dim.odims 0 2]);
    end
  end
end

% number of scans                                          
dim.S.sp = size(y.sp, length(dim.odims)+1);
dim.S.de = size(y.de, length(dim.odims)+1);

% number of echoes                                           
dim.E.sp = size(y.sp, length(dim.odims)+2);
dim.E.de = size(y.de, length(dim.odims)+2);

% default values 
arg.mask.est = [];
arg.mask.disp = [];
arg.mask.noise = [];

arg.nu.kap = [];
arg.nu.b0f = [];
arg.nu.b0s = [];
arg.nu.r2pf = [];
arg.nu.r2ps = [];

arg.wght.sp = ones(dim.S.sp, dim.E.sp);
arg.wght.de = ones(dim.S.de, dim.E.de);

arg.thresh = 0.05;

arg.x0.m0 = [];
arg.x0.ff = [];
arg.x0.t1f = [];
arg.x0.t1s = [];
arg.x0.t2f = [];
arg.x0.t2s = [];
arg.x0.kfs = [];

arg.meth.init = 'krr';
arg.meth.iter = 'pgpm';

arg.dist.x.m0.supp    = [];
arg.dist.x.m0.nsamp   = 1;
arg.dist.x.m0.prior   = 'unif';
arg.dist.x.ff.supp    = [-0.1 0.4].';
arg.dist.x.ff.nsamp   = 30;
arg.dist.x.ff.prior   = 'unif';
arg.dist.x.t1f.supp   = [50 700].';
arg.dist.x.t1f.nsamp  = 30;
arg.dist.x.t1f.prior  = 'logunif';
arg.dist.x.t1s.supp   = [700 2000].';
arg.dist.x.t1s.nsamp  = 30;
arg.dist.x.t1s.prior  = 'logunif';
arg.dist.x.t2f.supp   = [5 50].';
arg.dist.x.t2f.nsamp  = 30;
arg.dist.x.t2f.prior  = 'logunif';
arg.dist.x.t2s.supp   = [50 300].';
arg.dist.x.t2s.nsamp  = 30;
arg.dist.x.t2s.prior  = 'logunif';
arg.dist.x.kfs.supp   = [0.0001 0.05].';
arg.dist.x.kfs.nsamp  = 30;
arg.dist.x.kfs.prior  = 'logunif';

arg.dist.nu.kap.supp  = [0.5 2].';
arg.dist.nu.b0f.supp  = [0 0].';
arg.dist.nu.b0s.supp  = [0 0].';
arg.dist.nu.r2pf.supp = [0 0].';
arg.dist.nu.r2ps.supp = [0 0].';

arg.kmean.C = 10;
arg.kmean.opt = {...
  'EmptyAction', 'singleton',...
  'MaxIter', 1000};

arg.rff.snr = 0.1;
arg.rff.std = [];
arg.rff.len = [];
arg.rff.c = 2^0;
arg.rff.H = 10^4;
arg.rff.K = 10^6;

arg.train = [];

arg.inv.krr = 10^-8;
arg.inv.lm = 10^-8;

arg.precon.update = 5;
arg.precon.hessreg = 10^-10;

arg.line.step0 = 1;
arg.line.prop = 0.5;

arg.boxcon.m0 = [-Inf Inf].';
arg.boxcon.ff = [eps 0.40].';
arg.boxcon.t1f = [10 1000].';
arg.boxcon.t1s = [100 3000].';
arg.boxcon.t2f = [1 50].';
arg.boxcon.t2s = [10 300].';
arg.boxcon.kfs = [eps 0.1].';

arg.reg.m0.pot = {'hyper3', 2^-2};
arg.reg.m0.beta = [];
arg.reg.ff.pot = {'hyper3', 2^-6};
arg.reg.ff.beta = [];
arg.reg.t1f.pot = {'hyper3', 2^4};
arg.reg.t1f.beta = [];
arg.reg.t1s.pot = {'hyper3', 2^5};
arg.reg.t1s.beta = [];
arg.reg.t2f.pot = {'hyper3', 2^1};
arg.reg.t2f.beta = [];
arg.reg.t2s.pot = {'hyper3', 2^2};
arg.reg.t2s.beta = [];
arg.reg.kfs.pot = {'hyper3', 2^-7};
arg.reg.kfs.beta = [];

arg.stop.iter = 100;
arg.stop.wghtx = [];
arg.stop.tolx = 10^-7;

arg.bool.exchg  = 0;
arg.bool.mag.sp = 1;
arg.bool.mag.de = 1;
arg.bool.chat   = 0;
arg.bool.norm   = 1;
arg.bool.reset  = 1;
arg.bool.rfftst = 0;
arg.bool.nuclip = 1;
arg.bool.reg    = 0;
arg.bool.precon = 1;
arg.bool.disp   = 0;

arg.disp.m0 = [];
arg.disp.ff = [0 0.4];
arg.disp.t1f = [0 1000];
arg.disp.t1s = [0 2000];
arg.disp.t2f = [0 50];
arg.disp.t2s = [0 200];
arg.disp.kfs = [0 0.05];

% substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% if no display mask specified, use all voxels
if isempty(arg.mask.disp)
  arg.mask.disp = true(dim.odims);
  arg.mask.est  = true(dim.odims);
% else if no estimation mask specified, dilate mask.disp
elseif isempty(arg.mask.est)
  arg.mask.est = imdilate(arg.mask.disp, strel('disk', 10));
% else make sure that display mask tighter than estimation mask
elseif sum(arg.mask.disp & ~arg.mask.est) > 0
  error('Display mask has voxels outside estimation mask?!');
end

% noise mask
if isempty(arg.mask.noise)
  arg.mask.noise = ~arg.mask.est;
end

% set V to number of voxels in mask.est
dim.V = numel(arg.mask.est(arg.mask.est));

% check for spgr complex data
if ~isreal(y.sp) && arg.bool.mag.sp
  warn('\nDetected complex SPGR data but bool.mag.sp set true!? Reverting to default false.');
  arg.bool.mag.sp = false;
elseif ~isempty(y.sp) && isreal(y.sp) && ~arg.bool.mag.sp
  warn('\nDetected pure-real SPGR data; recommend setting bool.mag.sp true for speed.');
end

% check for dess complex data
if ~isreal(y.de) && arg.bool.mag.de
  warn('\nDetected complex DESS data but bool.mag.de set true!? Reverting to default false.');
  arg.bool.mag.de = false;
elseif ~isempty(y.de) && isreal(y.de) && ~arg.bool.mag.de
  warn('\nDetected pure-real DESS data; recommend setting bool.mag.de true for speed.');
end

% omit spgr scans with zero weight
keep.sp = arg.wght.sp(:,1)>0;
if ~all(keep.sp)
  y.sp = y.sp(:,:,keep.sp,:);
  P.sp.tr = P.sp.tr(keep.sp);
  P.sp.te = P.sp.te(keep.sp,:);
  P.sp.aex = P.sp.aex(keep.sp);
  arg.wght.sp = arg.wght.sp(keep.sp,:);
  dim.S.sp = size(arg.wght.sp,1);
  if arg.bool.chat
    warn('\nDetected zero-weight SPGR data: omitting %u datasets.',...
      sum(~keep.sp)*dim.E.sp);
  end
end

% omit dess scans with zero weight
keep.de = arg.wght.de(:,1)>0;
if ~all(keep.de)
  y.de = y.de(:,:,keep.de,:);
  P.de.tr = P.de.tr(keep.de);
  P.de.te = P.de.te(keep.de,:);
  P.de.aex = P.de.aex(keep.de);
  arg.wght.de = arg.wght.de(keep.de,:);
  dim.S.de = size(arg.wght.de,1);
  if arg.bool.chat
    warn('\nDetected zero-weight DESS data: omitting %u datasets.',...
      sum(~keep.de)*dim.E.de);
  end
end

% total number of datasets
dim.Dsp = dim.S.sp * dim.E.sp;
dim.Dde = dim.S.de * dim.E.de;
dim.D   = sum([dim.Dsp dim.Dde]);   

% vectorize data
y = cat(length(dim.odims)+1,...
  reshape(y.sp, [dim.odims dim.Dsp]),...
  reshape(y.de, [dim.odims dim.Dde]));                                  % [(odims) D]
if strcmp(arg.meth.init, 'krr')
  n = masker(y, arg.mask.noise);                                        % [V_bg D]
  n = col(transpose(n));                                                % [DV_bg]
end
y = masker(y, arg.mask.est);                                            % [V D]
y = col(transpose(y));                                                  % [DV]

% trick:  normalize data by median of sum of magnitude non-background values (over datasets)
%         effective regularization strength is then scale-invariant
if arg.bool.norm
  tmp = sum(reshape(abs(y), [dim.D dim.V]), 1);
  scale = median(tmp(tmp > arg.thresh * max(col(tmp))));
  y = div0(y,scale);
  if ~isempty(arg.x0.m0)
    arg.x0.m0 = div0(arg.x0.m0,scale);
  elseif strcmp(arg.meth.init, 'krr')
    n = div0(n,scale);
  end
end

% if nu fields unspecified, set to default values
if isempty(arg.nu.kap)
  arg.nu.kap = ones(dim.odims);
end
if isempty(arg.nu.b0f) || isempty(arg.nu.b0s) % either b0f/s empty
  if dim.S.de>0 && ~(arg.bool.mag.sp && arg.bool.mag.de) 
    warn('Using complex dess data w/o b0 map(s) will likely induce bias due to phase-mismatch.');
  end
  if ~isempty(arg.nu.b0f) % b0f but not b0s given
    warn('Fast- but not slow-compartment b0 given? Setting slow b0 to fast b0...');
    arg.nu.b0s = arg.nu.b0f;
  elseif ~isempty(arg.nu.b0s) % b0s but not b0f given
    warn('Slow- but not fast-compartment b0 given? Setting fast b0 to slow b0...');
    arg.nu.b0f = arg.nu.b0s;
  else % both empty
    arg.nu.b0f = zeros(dim.odims);
    arg.nu.b0s = zeros(dim.odims);
  end
end
if isempty(arg.nu.r2pf) || isempty(arg.nu.r2ps) % either r2pf/s empty
  if dim.S.de>0 && nnz(P.de.te(:,1)-P.de.te(:,2))>0
    warn('Using asymmetric dess echo times w/o r2p map(s) will induce bias.');
  end
  if ~isempty(arg.nu.r2pf) % r2pf but not r2ps given
    warn('Fast- but not slow-compartment r2p given? Setting slow to zero...');
    arg.nu.r2ps = zeros(dim.odims);
  elseif ~isempty(arg.nu.r2ps) % r2ps but not rp2f given
    warn('Slow- but not fast-compartment r2p given? Setting fast r2p to slow r2p...');
    arg.nu.r2pf = arg.nu.r2ps;
  else % both empty
    arg.nu.r2pf = zeros(dim.odims);
    arg.nu.r2ps = zeros(dim.odims);
  end
end

% vectorize nu fields
arg.nu.kap  = arg.nu.kap(arg.mask.est);                                 % [V]
arg.nu.b0f  = arg.nu.b0f(arg.mask.est);                                 % [V]
arg.nu.b0s  = arg.nu.b0s(arg.mask.est);                                 % [V]  
arg.nu.r2pf = arg.nu.r2pf(arg.mask.est);                                % [V]
arg.nu.r2ps = arg.nu.r2ps(arg.mask.est);                                % [V]

% vectorize data weights
w = [...
  col(arg.wght.sp);...
  col(arg.wght.de)];                                                    % [D]

% if neglecting exchange, remove corresponding fields
if ~arg.bool.exchg
  arg.x0 = rmfield(arg.x0, 'kfs');
  arg.dist.x = rmfield(arg.dist.x, 'kfs');
  arg.boxcon = rmfield(arg.boxcon, 'kfs');
  arg.reg = rmfield(arg.reg, 'kfs');
  arg.disp = rmfield(arg.disp, 'kfs');
end

% warn user if attempt is made to constrain m0 estimation
if ~all(isinf(arg.boxcon.m0))
  warn('Unconstrained m0 estimation: ignoring arg.boxcon.m0.* settings.');
  arg.boxcon.m0 = [-Inf Inf].';
end

% save fieldnames
field.x = fieldnames(arg.x0);
field.nu = fieldnames(arg.nu);

% to enumerate conveniently, switch from struct to cell
arg.nu = struct2cell(arg.nu);                                           % {N} cell of [V] arrays
arg.dist.x = struct2cell(arg.dist.x);                                   % {L} cell of [1] structs
arg.dist.nu = struct2cell(arg.dist.nu);                                 % {N} cell of [1] structs
arg.x0 = struct2cell(arg.x0);                                           % {L} cell of [(odims)] arrays
arg.boxcon = struct2cell(arg.boxcon);                                   % {L} cell of [2] arrays
arg.disp = struct2cell(arg.disp);                                       % {L} cell of [1 2] arrays
dim.L = length(arg.x0);
dim.N = length(arg.nu);

% if any x0 fields left unspecified, initialization necessary
tmp = 0;
for l = 1:dim.L
  tmp = tmp + ~isempty(arg.x0{l});
end
if tmp<dim.L
  % choose initialization method
  switch arg.meth.init
    case 'vpm'
      % warn user if attempt is made to constrain m0 ml estimation
      if ~isempty(arg.dist.x{1}.supp) || arg.dist.x{1}.nsamp~=1
        warn('\nVarPro ML estimation fixes unity m0: ignoring arg.dist.x.m0.* settings.');
        arg.dist.x{1}.supp  = [];
        arg.dist.x{1}.nsamp   = 1;
        arg.dist.x{1}.prior   = 'unif';
      end
      arg.dist.x{1}.supp = [1 1]';
      
      % use only one k-means cluster if all nu and preset x0 fields uniform
      tmp = true;
      for n = 1:dim.N
        if min(arg.nu{n})~=max(arg.nu{n})
          tmp = false;
          break;
        end
      end
      if tmp
        for l = 1:dim.L
          if ~isempty(arg.x0{l})
            tmp2 = minmax(masker(arg.x0{l}, arg.mask.est));
            if min(tmp2)~=max(tmp2)
              tmp = false;
              break;
            end
          end
        end
      end
      if tmp
        fprintf('Detected uniform known maps: setting kmean.C to 1.\n');
        arg.kmean.C = 1;
      end

      % max-likelihood estimation via variable projection method
      if arg.bool.chat
        fprintf('\n========================================================================');
        fprintf('\nInitialization via VarPro...');
        fprintf('\n========================================================================\n');
      end
      [arg.x0, t.init] = varpro(...
        y, arg.x0, arg.nu, P, w,...
        arg.kmean, arg.dist.x, dim, arg.bool, arg.mask.est);            % {L} cell w/ [V] blocks
    case 'krr'
      % check max sig for unity m0
      if arg.rff.snr>1
        error('max unity-m0 signal cannot exceed 1!');
      elseif arg.rff.snr<=eps
        error('max unity-m0 signal must be positive and less than 1!');
      elseif arg.rff.snr>0.3
        tmp = strcat(...
          'Detected max unity-m0 signal >0.3: are you sure?\n',...
          'Overestimating may induce improper m0 sampling and krr errors...');
        warn(tmp);
      elseif arg.rff.snr<0.01
        tmp = strcat(...
          'Detected max unity-m0 signal <0.01: are you sure?\n',...
          'Underestimating will cause inefficient m0 sampling and krr underperformance...');
        warn(tmp);
      end
      
      % set training data noise std dev
      if isempty(arg.rff.std)
        if isempty(n)
          if arg.bool.norm
            arg.rff.std = div0(3.8607e-4,scale);
          else
            arg.rff.std = 3.8607e-4;
          end
          warn('No noise mask given: rff.std set automatically to %0.8f.', arg.rff.std);
        else
          if isreal(n)
            % n is assumed ~rayleigh(\sigma)
            % rff.std is a (slightly biased) estimate of \sigma
            tmp = sum(n.^2);
            tmp = div0(tmp,2*length(n));
            arg.rff.std = sqrt(tmp);
          else
            error('todo: estimate rff.std for complex data');
          end
        end
      end
      
      % set kernel length scales 
      if isempty(arg.rff.len)
        tmp = reshape(y, [dim.D dim.V]);
        tmp = [tmp; transpose([arg.nu{:}])];
        arg.rff.len = mean(tmp,2);
        arg.rff.len = max(arg.rff.len,eps);
      end

      % kernel regression
      if arg.bool.chat
        fprintf('========================================================================');
        fprintf('\nInitialization via Kernel Regression...');
        fprintf('\n========================================================================\n');
      end
      [arg.x0, t.init] = krr(...
        y, arg.x0, arg.nu, P, w,...
        arg.rff, arg.train, arg.dist, arg.inv.krr,...
        dim, arg.bool, arg.mask.est);                                 % {L} cell w/ [V] blocks
    otherwise
      error('Unknown initialization method requested!');
  end
  if arg.bool.chat  
    fprintf('========================================================================');
    fprintf('\n                                    ...done in %0.3f seconds.', t.init);
    fprintf('\n========================================================================\n');
  end
else
  for l = 1:dim.L
    arg.x0{l} = masker(arg.x0{l}, arg.mask.est);
  end
end

if arg.bool.reg
  % if not preset, set regularizer strengths proportional to D
  if isempty(arg.reg.m0.beta)
    arg.reg.m0.beta = dim.D * 2^-26;
  end
  if isempty(arg.reg.ff.beta)
    arg.reg.ff.beta = dim.D * 2^-26;
  end
  if isempty(arg.reg.t1f.beta)
    arg.reg.t1f.beta = dim.D * 2^-22;
  end
  if isempty(arg.reg.t1s.beta)
    arg.reg.t1s.beta = dim.D * 2^-21;
  end
  if isempty(arg.reg.t2f.beta)
    arg.reg.t2f.beta = dim.D * 2^-25;
  end
  if isempty(arg.reg.t2s.beta)
    arg.reg.t2s.beta = dim.D * 2^-23;
  end
  if dim.L>6 && isempty(arg.reg.kfs.beta)
    arg.reg.kfs.beta = dim.D * 2^-28;
  end

  % to enumerate conveniently, switch from struct to cell
  arg.reg = struct2cell(arg.reg);                                         % {L}

  % regularizer objects
  for l = 1:dim.L
    arg.reg{l}.R = Reg1(arg.mask.est,...
      'pot_arg', arg.reg{l}.pot,...
      'beta', arg.reg{l}.beta,...
      'type_penal', 'mat');
  end
end

% set iterative stopping criterion weights
arg.stop.wghtx = NaN(dim.L,1);
for l=1:dim.L
  if l==1
    tmp = arg.x0{l};
    tmp = median(tmp(tmp > arg.thresh * max(col(tmp))));
  else
    tmp = mean(arg.boxcon{l});
  end
  arg.stop.wghtx(l) = div0(1,tmp);
end

% choose iterative method
switch arg.meth.iter
  case 'plm'
    if arg.bool.chat
      fprintf('\n========================================================================');
      fprintf('\nIterative estimation via projected levenberg-marquardt...');
      fprintf('\n========================================================================\n');
    end
    [tmp, t.iter] = plm(...
      y, arg.x0, arg.nu, P, w, arg.inv.lm, arg.boxcon,...
      arg.reg, arg.stop, dim, arg.bool, arg.mask.est, arg.disp);        % {L} cell w/ [V] blocks
  case 'pgpm'
    if arg.bool.chat
      fprintf('\n========================================================================');
      fprintf('\nIterative estimation via gradient projection method...');
      fprintf('\n========================================================================\n');
    end
    [tmp, t.iter] = pgpm(...
      y, arg.x0, arg.nu, P, w, arg.kmean, arg.precon,...
      arg.line, arg.boxcon, arg.reg, arg.stop, dim, arg.bool);          % {L} cell w/ [V] blocks
  otherwise
    error('Unknown iterative method!');
end
if arg.bool.chat
  fprintf('\n========================================================================');
  fprintf('\n                        ...done in %0.3f seconds.', t.iter);
  fprintf('\n========================================================================\n\n');
end

% trick: rescale m0 for output
if arg.bool.norm
  arg.x0{1} = arg.x0{1} * scale;
  tmp{1} = tmp{1} * scale;
end

% postprocessing
for i = 1:length(field.x)
  % embed for output
  x.init.(field.x{i}) = embed(arg.x0{i}, arg.mask.est);                 % [(odims)]
  x.iter.(field.x{i}) = embed(tmp{i}, arg.mask.est);                    % [(odims)]
  
  % apply mask.disp for display
  x.init.(field.x{i})(~arg.mask.disp) = 0;
  x.iter.(field.x{i})(~arg.mask.disp) = 0;
end
end
