  function [x, t, t2] = mri_mese_mwf(y, P, varargin)
%|function [x, t, t2] = mri_mese_mwf(y, P, varargin)
%|
%|  myelin water fraction mapping from multi-echo spin echo data
%|    estimates t2 distribution using nonnegative least-squares
%|    accomodates incomplete refocusing via extended phase graph formalism
%|
%|  inputs
%|    y         [(dim.o) dim.e] coil-combined magnitude image data
%|    P         [1x1 struct]    scan parameters
%|     .ex                        excitation parameters
%|        .a    [1]                 nominal flip angle                                        rad 
%|        .ph   [1]                 phase relative to +x axis                                 rad
%|     .ref                       refocusing parameters 
%|         .a   [dim.e]             refocusing nominal flip angle                             rad
%|         .ph  [dim.e]             refocusing phase relative to +x axis                      rad
%|     .ncyc    [dim.e 2]         across-voxel phase cycles
%|                                  imparted by gradient crusher pairs        
%|     .tr      [1]               repetition times                                            ms
%|     .te      [dim.e]           echo times                                                  ms
%|
%|  options
%|    mask      [1x1 struct]    binary object masks
%|     .disp    [(dim.o)]         over which to display images          def: true(dim.o)
%|     .est     [(dim.o)]         over which to perform estimation      def: imdilate(mask.disp)
%|    nu        [1x1 struct]    known parameters
%|     .t1      [(dim.o)]         spin-lattice relaxation time map      def: 1000*ones(dim.o) ms
%|     .kap     [(dim.o)]         flip-angle scale map                  def: ones(dim.o)
%|    wght      [dim.e]         dataset weights                         def: ones(dim.e,1)
%|    thresh    [1]             frac of max sig deciding 'background'   def: 0.05
%|    dist      [1x1 struct]    t2 sampling distribution object
%|     .boxcon  [2]               [lb ub] distribution endpoints        def: [10 3000]        ms
%|     .nsamp   [1]               number of dictionary samples          def: 100 
%|     .prior   {1}               ('unif', 'logunif') distribution      def: 'logunif'
%|    reg       [1]             l2 regularization parameter             def: 2^-30
%|    kmean     [1x1 struct]    kmeans object to pool nu
%|     .C       [1]               number of known clusters              def: 100
%|     .opt     {1 2*?}           optional arguments to kmeans          def: see below
%|    nnls      {1 2*?}         optional arguments to lsqnonneg         def: see below
%|    comp      [1x1 struct]    compartmental discretization object
%|     .name    {dim.comp}        water compartment fraction names      def: {'mw','iew','fw'}
%|     .boxcon  {dim.comp}        [lb ub] compartment t2 endpoints      def: see below        ms
%|    bool      [1x1 struct]    boolean variables
%|     .chat    false|true        verbosity                             def: true
%|     .scale   false|true        normalize data for scale-invar reg    def: true
%|     .clust   false|true        run nnls over clustered voxels        def: true
%|     .reg     false|true        use l2 regularization                 def: true
%|     .pool    false|true        use parallel workers                  def: false
%|     .norm    false|true        normalize w.r.t. total water content  def: true
%|
%|  outputs
%|    x         [1x1 struct]    compartmental fraction estimates
%|     .tot     [(dim.o)]         total spin density estimate
%|     .*       [(dim.o)]         fieldnames(x) given by comp.name
%|    t         [1]             mese nnls run time                                            s
%|    t2        [1x1 struct]    t2 distribution object
%|     .samp    [dim.s]           t2 sample values                                            ms
%|     .dist    [(dim.o) dim.s]   unit-norm t2 distribution estimates
%|  
%|  copyright 2017-8, gopal nataraj, university of michigan
%|
%|  version control
%|    2017-11-02                original
%|    2017-11-05                clustering now optional
%|    2018-02-22                optional l2 regularization; return t2 distribution

% dimensions
tmp = size(y);
dim.o = tmp(1:end-1);
dim.e = tmp(end);

% default values 
arg.mask.disp = [];
arg.mask.est = [];
arg.nu.t1 = [];
arg.nu.kap = [];
arg.wght = ones(dim.e,1);
arg.thresh = 0.05;
arg.dist.boxcon = [10 3000].';
arg.dist.nsamp = 100;
arg.dist.prior = 'logunif';
arg.reg = 2^(-30);
arg.kmean.C = 100;
arg.kmean.opt = {...
  'EmptyAction', 'singleton',...
  'MaxIter', 1000};
arg.nnls = {};
arg.comp.name = {'mw'; 'iew'; 'fw'};
arg.comp.boxcon = {...
  [15 40].';...
  [40 200].';...
  [200 Inf].'
};
arg.bool.chat = 1;
arg.bool.scale = 1;
arg.bool.clust = 1;
arg.bool.reg = 1;
arg.bool.pool = 0;
arg.bool.norm = 1;

% substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% throw error if any compartment is named 'tot'
for i = 1:length(arg.comp.name)
  if strcmp(arg.comp.name{i}, 'tot')
    error('Compartment name ''tot'' reserved for total spin density.');
  end
end

% if no display mask specified, use all voxels
if isempty(arg.mask.disp)
  arg.mask.disp = true(dim.o);
  arg.mask.est = true(dim.o);
% else if no estimation mask specified, dilate mask.disp
elseif isempty(arg.mask.est)
  arg.mask.est = imdilate(arg.mask.disp, strel('disk', 10));
% else ensure that display mask tighter than estimation mask
elseif sum(arg.mask.disp & ~arg.mask.est) > 0
  error('Display mask has voxels outside estimation mask?!');
end

% set dim.v to number of voxels in mask.est
dim.v = numel(arg.mask.est(arg.mask.est));

% if input images complex, take magnitude and warn user
if ~isreal(y)
  y = abs(y);
  if arg.bool.chat
    warn('Detected complex y but lsqnonneg requires pure-real inputs: using abs(y).');
  end
end

% vectorize data
y = masker(y, arg.mask.est);                                            % [dim.v dim.e]

% omit scans with zero weight
tmp = arg.wght>0;
if ~all(tmp)
  y = y(:,tmp);
  P.ref.a = P.ref.a(tmp);
  P.ref.ph = P.ref.ph(tmp);
  P.ncyc = P.ncyc(tmp,:);
  arg.wght = arg.wght(tmp);
  if arg.bool.chat
    warn('Detected zero-weight data: omitting %u datasets.', sum(~tmp));
  end
end

% scale reg parameter by square of median of first-echo magnitude non-background values
% effective regularization strength is then roughly scale-invariant
if arg.bool.scale && arg.bool.reg
  tmp = y(:,1);
  scale = median(tmp(tmp > arg.thresh*max(tmp)));
  arg.reg = arg.reg * scale^2;
end

% if nu fields unspecified, set to defaults and warn user
if isempty(arg.nu.t1)
  arg.nu.t1 = 1000*ones(dim.o);
  if arg.bool.chat
    warn('nu.t1 not specified: setting to 1000 everywhere.');
  end
end
if isempty(arg.nu.kap)
  arg.nu.kap = ones(dim.o);
  if arg.bool.chat
    warn('nu.kap not specified: setting to 1 everywhere.');
  end
end

% vectorize nu fields
arg.nu.t1 = arg.nu.t1(arg.mask.est);                                    % [dim.v]
arg.nu.kap = arg.nu.kap(arg.mask.est);                                  % [dim.v]

% use only one k-means cluster if all nu fields uniform
tmp = fieldnames(arg.nu);
tmp2 = true;
for n = 1:length(tmp)
  if min(arg.nu.(tmp{n}))~=max(arg.nu.(tmp{n}))
    tmp2 = false;
    break;
  end
end
if tmp2
  arg.kmean.C = 1;
  if arg.bool.chat
    fprintf('Detected uniform known maps: setting kmean.C to 1.\n');
  end
end

% set t2 distribution sample locations
switch arg.dist.prior
  case 'unif'
    t2.samp = col(linspace(...
      arg.dist.boxcon(1), arg.dist.boxcon(2), arg.dist.nsamp));
  case 'logunif'
    t2.samp = col(logspace(...
      log10(arg.dist.boxcon(1)), log10(arg.dist.boxcon(2)), arg.dist.nsamp));
  otherwise
    error('Unknown t2 distribution prior option.');
end
dim.s = arg.dist.nsamp;

% nonnegative least-squares estimation of t2 distribution
if arg.bool.chat
  fprintf('\n==================================================');
  fprintf('\nnnls estimation of t2 distribution...');
  fprintf('\n==================================================\n');
end
[t2.dist, t] = mese_nnls(...
  y, P, arg.nu, arg.wght, t2.samp,...
  arg.reg, arg.kmean, arg.nnls, arg.bool, dim);                         % [dim.v dim.s]
if arg.bool.chat  
  fprintf('=============================================================');
  fprintf('\n                                    ...done in %0.3f seconds.', t);
  fprintf('\n=============================================================\n');
end

% total spin density
x.tot = sum(t2.dist,2);                                                 % [dim.v]

% pool into compartments
tmp = col(arg.comp.name);
for i = 1:length(tmp)
  tmp2 = arg.comp.boxcon{i}(1)<t2.samp;
  tmp2 = tmp2 & t2.samp<=arg.comp.boxcon{i}(2);                         % [dim.s]
  x.(tmp{i}) = sum(t2.dist(:,tmp2),2);                                  % [dim.v]
  if arg.bool.norm
    x.(tmp{i}) = div0(x.(tmp{i}), x.tot);
  end
end

% normalize and embed distribution
if nargout>2
  t2.dist = bsxfun(@rdivide, t2.dist, x.tot);                           % [dim.v dim.s]
  t2.dist = embed(t2.dist, arg.mask.est);                               % [(dim.o) dim.s]
end

% embed for output
tmp = fieldnames(x);
for i = 1:length(tmp)
  x.(tmp{i}) = embed(x.(tmp{i}), arg.mask.est);                         % [(dim.o)]
  x.(tmp{i})(~arg.mask.disp) = 0;
end
end
  