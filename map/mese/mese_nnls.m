  function [t2dist, t] = mese_nnls(y, P, nu, w, t2, reg, kmean, nnls, bool, dim)
%|function [t2dist, t] = mese_nnls(y, P, nu, w, t2, reg, kmean, nnls, bool, dim)
%|
%|  nonnegative least-squares estimation of t2 distribution
%|    preclusters known parameter maps and computes separate system matrices for each cluster.
%|    
%|  inputs
%|    y         [dim.v dim.e]   coil-combined magnitude image data
%|    P         [1x1 struct]    scan parameters
%|     .ex                        excitation parameters
%|        .a    [1]                 nominal flip angle                                        rad 
%|        .ph   [1]                 phase relative to +x axis                                 rad
%|     .ref                       refocusing parameters 
%|         .a   [dim.e]             refocusing nominal flip angle                             rad
%|         .ph  [dim.e]             refocusing phase relative to +x axis                      rad
%|    nu        [1x1 struct]    known parameters
%|     .t1      [dim.v]           spin-lattice relaxation time map      
%|     .kap     [dim.v]           flip-angle scale map                
%|    w         [dim.e]         dataset weights
%|    t2        [dim.s]         t2 values at which to estimate t2 dist                        ms
%|    reg       [1]             l2 regularization parameter         
%|    kmean     [1x1 struct]    kmeans object to pool nu
%|     .C       [1]               number of known clusters              
%|     .opt     {1 2*?}           optional arguments to kmeans      
%|    nnls      {1 2*?}         optional arguments to lsqnonneg     
%|    bool      [1x1 struct]    boolean variables
%|     .chat    false|true        verbosity                       
%|     .clust   false|true        run nnls over clustered voxels      
%|     .pool    false|true        use parallel workers              
%|     .norm    false|true        normalize w.r.t. total water content 
%|    dim       [1x1 struct]    object containing dimension info
%|
%|  outputs
%|    t2dist    [dim.v dim.s]   unnormalized t2 distribution estimate 
%|                                sum(t2dist,2) gives spin density est
%|    t         [1]             run time                                                      s
%|
%|  copyright 2017, gopal nataraj, university of michigan
%|
%|  version control
%|    2017-11-02                original
%|    2017-11-05                clustering now optional

% start clock
tic;

% compute system matrix for unity spin density
m0 = ones(dim.s,1);

% instantiate t2 distribution estimate
t2dist = nan(dim.v, dim.s);

% data weighting matrix
sqrtW = spdiags(sqrt(w), 0, dim.e, dim.e);                              % [dim.e dim.e sparse]

% l2 regularization constraints
if bool.reg
  A_reg = diag(ones(dim.s,1)*sqrt(reg));                                % [dim.s dim.s]
end

% bool.clust==1: faster, but less accurate
if bool.clust
  if bool.chat
    fprintf('mese_nnls: computing separate dictionary per voxel cluster...\n');
  end
  
  % cluster nu maps
  [cl.idx, cl.cent] = kmeans([nu.t1 nu.kap], kmean.C, kmean.opt{:});    % [dim.v], [kmean.C 2]
  
  % for each cluster, construct system matrix and run lsqnonneg on voxels assigned to that cluster
  for c = 1:kmean.C
    if bool.chat
      fprintf('mese_nnls: cluster %2u of %2u...', c, kmean.C);
    end

    % system matrix fixed for cth cluster
    A = transpose(mese(...
      m0, ones(dim.s,1)*cl.cent(c,1), t2, ones(dim.s,1)*cl.cent(c,2),...
      P.ex, P.ref, P.ncyc, P.tr, P.te, 'bool.mag', true));              % [dim.e dim.s]

    % temporary variables to minimize memory requirement of parallel workers
    yc = transpose(y(cl.idx==c,:));                                     % [dim.e dim.vc]
    t2distc = nan(dim.s,size(yc,2));                                    % [dim.s dim.vc]

    % apply data weighting
    A = sqrtW * A;
    yc = sqrtW * yc;
    
    % impose l2 regularization constraints
    if bool.reg
      A = [A; A_reg];                                                   % [dim.e+dim.s dim.s]
      yc = [yc; zeros(dim.s,size(yc,2))];                               % [dim.e+dim.s dim.vc]
    end

    % bool.pool==1: parallelize across voxels within cth cluster
    if bool.pool 
      parfor i = 1:sum(cl.idx==c)
        t2distc(:,i) = lsqnonneg(A, yc(:,i), nnls{:});
      end
    else
      for i = 1:sum(cl.idx==c)
        t2distc(:,i) = lsqnonneg(A, yc(:,i), nnls{:});
      end
    end

    % store cth cluster results
    t2dist(cl.idx==c,:) = transpose(t2distc);                           % [dim.vc dim.s]

    if bool.chat
      fprintf('mapped %6u of %6u voxels.\n', sum(cl.idx<=c), dim.v);
    end
  end
% bool.clust==0: slower, but more accurate
else
  if bool.chat
    fprintf('mese_nnls: computing separate dictionary per voxel...\n');
  end
  
  % store for speed
  if bool.reg
    y_reg = zeros(dim.s,1);
  end
  
  % bool.pool==1: parallelize across voxels 
  if bool.pool
    % temporary variables to minimize memory requirement of parallel workers
    tmp1 = ones(dim.s,1)*transpose(nu.t1);                                % [dim.s dim.v]
    tmp2 = ones(dim.s,1)*transpose(nu.kap);                               % [dim.s dim.v]
    tmp3 = bool.reg;
    
    parfor i = 1:dim.v
      % system matrix
      Ai = transpose(mese(...
        m0, tmp1(:,i), t2, tmp2(:,i),...
        P.ex, P.ref, P.ncyc, P.tr, P.te, 'mag', true));                   % [dim.e dim.s]
      
      % apply data weighting
      Ai = sqrtW * Ai;
      yi = sqrtW * transpose(y(i,:));                                     % [dim.e]
      
      % impose l2 regularization constraints
      if tmp3
        Ai = [Ai; A_reg];                                                 % [dim.e+dim.s dim.s]
        yi = [yi; y_reg];                                                 % [dim.e+dim.s]
      end
      
      % nonnegative least-squares
      tmp = lsqnonneg(Ai, yi, nnls{:});
      t2dist(i,:) = transpose(tmp);
    end
  else
    for i = 1:dim.v
      % system matrix
      Ai = transpose(mese(...
        m0, ones(dim.s,1)*nu.t1(i), t2, ones(dim.s,1)*nu.kap(i),...
        P.ex, P.ref, P.ncyc, P.tr, P.te, 'mag', true));                   % [dim.e dim.s]
      
      % apply data weighting
      Ai = sqrtW * Ai;
      yi = sqrtW * transpose(y(i,:));                                     % [dim.e]
      
      % impose l2 regularization constraints
      if bool.reg
        Ai = [Ai; A_reg];                                                 % [dim.e+dim.s dim.s]
        yi = [yi; y_reg];                                                 % [dim.e+dim.s]
      end
      
      % nonnegative least-squares
      tmp = lsqnonneg(Ai, yi, nnls{:});
      t2dist(i,:) = transpose(tmp);
      
      % optional: updates
      if bool.chat && rem(i,1000)==0
        fprintf('mese_nnls: mapped %u of %u voxels.\n', i, dim.v);
      end
    end
  end
end
t = toc;
end