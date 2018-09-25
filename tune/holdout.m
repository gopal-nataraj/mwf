% script holdout.m
% hyperparameter tuning for krr option of mri_multicomp_map.m
%
% copyright 2018, gopal nataraj, university of michigan
%
% version control
%   2018-02-23      adapted from m0,t1,t2 estimation version

% irt 
if (~exist('irtdir', 'var'))
  curdir = pwd;
  cd ../../../../../irt;
  irtdir = pwd;
  setup();
  cd(curdir);
end

% mapping
addpath('../../model/spgr');
addpath('../../model/dess');
addpath('../../map');
addpath('../../etc');

% noise standard deviation
c.sig = 3.8607e-4; 

% latent parameter distribution
dist.x.m0.supp      = [0.001 1];
dist.x.m0.prior     = 'unif';
dist.x.ff.supp      = [0.03 0.21].';  % needs to be positive
dist.x.ff.prior     = 'unif';
dist.x.t1f.supp     = [50 700].';
dist.x.t1f.prior    = 'logunif';
dist.x.t1s.supp     = [700 2000].';
dist.x.t1s.prior    = 'logunif';
dist.x.t2f.supp     = [5 50].';
dist.x.t2f.prior    = 'logunif';
dist.x.t2s.supp     = [50 300].';
dist.x.t2s.prior    = 'logunif';

% known parameter distribution
dist.nu.kap.supp    = [0.5 2];
dist.nu.kap.prior   = 'lognormal';    % log(nu) is normal
dist.nu.kap.logmu   = 0;              % log(nu) mean
dist.nu.kap.logsd   = 0.2;            % log(nu) std dev
dist.nu.b0f.supp    = [0 0].';
dist.nu.b0f.prior   = 'unif';
dist.nu.b0s.supp    = [0 0].';
dist.nu.b0s.prior   = 'unif';
dist.nu.r2pf.supp   = [0 0].';
dist.nu.r2pf.prior  = 'unif';
dist.nu.r2ps.supp   = [0 0].';
dist.nu.r2ps.prior  = 'unif';

% header options
header.tru = 0;                       % initialize with truth (to test)
header.sv = 1;                        % save holdout results
header.im = 0;                        % show holdout plots
header.pr = 0;                        % print holdout plots

% holdout options
hld.lam = 2.^col(linspace(-1.5,1.5,31));
hld.rho = 2.^col(linspace(-35,-5,31));
hld.V = 10^5;

% sample latent parameters
field = fieldnames(dist.x);
for l = 1:length(field)
  tmp = dist.x.(field{l}).supp;
  switch dist.x.(field{l}).prior
    case 'unif'
      x.(field{l}) = random('unif', tmp(1), tmp(2), [hld.V 1]);
      
    case 'logunif'
      tmp = log(tmp);
      tmp2 = random('unif', tmp(1), tmp(2), [hld.V 1]);
      x.(field{l}) = exp(tmp2);
      
    otherwise
      error('unknown dist for latent parameter %u.', l);
  end
  x.(field{l}) = reshape(x.(field{l}), 100, []);
end

% sample known parameters
field = fieldnames(dist.nu);
for k = 1:length(field)
  tmp = dist.nu.(field{k}).supp;
  switch dist.nu.(field{k}).prior
    case 'unif'
      nu.(field{k}) = random('unif', tmp(1), tmp(2), [hld.V 1]);
      
    case 'logunif'
      tmp = log(tmp);
      tmp2 = random('unif', tmp(1), tmp(2), [hld.V 1]);
      nu.(field{k}) = exp(tmp2);
      
    case 'lognormal'
      tmp2 = makedist('normal',...
        'mu', dist.nu.(field{k}).logmu,...
        'sigma', dist.nu.(field{k}).logsd);
      tmp2 = truncate(tmp2, log(tmp(1)), log(tmp(2)));
      tmp2 = random(tmp2, [hld.V 1]);
      nu.(field{k}) = exp(tmp2);
      
    otherwise
      error('unknown dist for known parameter %u.', k);
  end
  nu.(field{k}) = reshape(nu.(field{k}), 100, []);
end

% dess scan parameters
P.de.tr = [17.500; 30.197; 60.303];       % ms
P.de.te = ones(3,1) * [5.29 5.29];        % ms
P.de.aex = [32.995; 18.303; 15.147];      % deg
P.de.aex = P.de.aex * (pi/180);           % rad
S.de = size(P.de.te,1);
E.de = size(P.de.te,2);

% recon options
meth.init = 'krr';
meth.iter = 'pgpm';
rff.snr = 0.1;
rff.std = c.sig;
rff.len = [];
rff.c = [];
rff.H = 10^3;
rff.K = 10^6;
inv.krr = [];
stop.iter = 0;
bool.exchg = 0;
bool.mag.sp = 1;
bool.mag.de = 1;
bool.chat = 2;
bool.norm = 1;
bool.reset = 1;
bool.rfftst = 0;
bool.nuclip = 1;
bool.reg = 0;
bool.precon = 1;
bool.disp = 0;

% acquire noiseless dess data
y.de = NaN([size(x.m0), S.de, E.de]);
for s = 1:S.de
  [y.de(:,:,s,1), y.de(:,:,s,2)] = dess_2comp(...
    x.m0, x.ff, x.t1f, x.t1s, x.t2f, x.t2s,...
    nu.kap, nu.b0f, nu.b0s,...
    P.de.aex(s), P.de.tr(s), P.de.te(s,1), P.de.te(s,2),...
    'R2p_f', nu.r2pf,...
    'R2p_s', nu.r2ps,...
    'exchg', bool.exchg,...
    'mag', 0);
end

% add complex gaussian noise
n.de = c.sig * (randn(size(y.de)) + 1i*randn(size(y.de)));
ycc.de = y.de + n.de;

% use magnitude data as appropriate
if bool.mag.de, ycc.de = abs(ycc.de); end

% % set dist.x.m0.supp to match map setting
% dist.x.m0.supp = [eps div0(max(col(ycc.de)),rff.snr)].';

% initial estimate
if header.tru
  x0 = x;
else
  x0.m0 = [];
  x0.ff = [];
  x0.t1f = [];
  x0.t1s = [];
  x0.t2f = [];
  x0.t2s = [];
  x0.kfs = [];
end

% krr hyperparameter optimization
try
  load('holdout.mat');
catch
  % weighted nrmse measures
  nrmse.m0  = nan(length(hld.lam), length(hld.rho));
  nrmse.ff  = nan(length(hld.lam), length(hld.rho));
  nrmse.t1f  = nan(length(hld.lam), length(hld.rho));
  nrmse.t1s  = nan(length(hld.lam), length(hld.rho));
  nrmse.t2f  = nan(length(hld.lam), length(hld.rho));
  nrmse.t2s  = nan(length(hld.lam), length(hld.rho));
  nrmse.all = nan(length(hld.lam), length(hld.rho));  
  
  % parameter weighting
  w.m0  = [1 0 0 0 0 0].';
  w.ff  = [0 1 0 0 0 0].';
  w.t1f = [0 0 1 0 0 0].';
  w.t1s = [0 0 0 1 0 0].';
  w.t2f = [0 0 0 0 1 0].';
  w.t2s = [0 0 0 0 0 1].';
  w.all = ones(6,1)*(1/6);

  % record total holdout time
  t = 0;
  for h1 = 1:length(hld.lam)
    rff.c = hld.lam(h1);
    
    for h2 = 1:length(hld.rho)
      inv.krr = hld.rho(h2);

      % parameter estimation
      fprintf('\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
      fprintf(...
        'Holdout: (rff.c, inv.krr) = (2^%.2f, 2^%.2f)...\n',...
        log2(rff.c), log2(inv.krr));
      fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
      opt.map = {...
        'nu', nu,...
        'x0', x0,...
        'dist.x.ff', dist.x.ff,...
        'dist.x.t1f', dist.x.t1f,...
        'dist.x.t1s', dist.x.t1s,...
        'dist.x.t2f', dist.x.t2f,...
        'dist.x.t2s', dist.x.t2s,...
        'meth.init', meth.init,...
        'meth.iter', meth.iter,...
        'rff', rff,...
        'inv.krr', inv.krr,...
        'stop.iter', stop.iter,...
        'bool', bool
      };
      [xhat, time] = mri_multicomp_map(ycc, P, opt.map{:});
      t = t + time.init;
      
      % holdout nrmse of krr
      tmp = fieldnames(nrmse);
      for i = 1:length(tmp)
        opt.wnrmse = {...
          'wght', w.(tmp{i})
        };
        nrmse.(tmp{i})(h1,h2) = wnrmse(xhat.init, x, opt.wnrmse{:});
      end
    end
  end
  
  % save wmse
  if header.sv
    save('holdout.mat', 'nrmse', 'w', 'hld', 't');
  end
end

% find minima
tmp = fieldnames(nrmse);
for i = 1:length(tmp)
	[tmp2, idx.lam{i}] = min(nrmse.(tmp{i}),[],1);
  [tmp3, idx.rho{i}] = min(tmp2);
  fprintf('\n%3s-weighted nrmse minimized at (lam,rho) = (2^%.2f,2^%.2f) with value %0.4f.\n',...
    tmp{i},...
    log2(hld.lam(idx.lam{i}(idx.rho{i}))),...
    log2(hld.rho(idx.rho{i})),...
    nrmse.(tmp{i})(idx.lam{i}(idx.rho{i}), idx.rho{i}));
end

% weighted nrmse plots
if header.im
  % tick labels
  tick.x = log2(hld.lam(1:5:end));
  tick.y = log2(hld.rho(1:5:end));
  for i = 1:length(tick.x)
    ticklab.x{i} = sprintf('$2^{%0.1f}$', tick.x(i));
  end
  for i = 1:length(tick.y)
    ticklab.y{i} = sprintf('$2^{%d}$', tick.y(i));
  end
  
  % dynamic ranges
  tmp = fieldnames(nrmse);
  for i = 1:length(tmp)
    tmp2 = minmax(nrmse.(tmp{i}));
    tmp2(1) = floor(100*tmp2(1))/100;
    tmp2(2) = ceil(100*tmp2(2))/100;
    dyn.(tmp{i}) = tmp2.';
  end
  
  % plots
  tmp = fieldnames(nrmse);
  for i = 1:length(tmp)
    figure;
    hold on;
    im(log2(hld.lam), log2(hld.rho), nrmse.(tmp{i}), dyn.(tmp{i}), 'cbar', ' ');
    lim.x = xlim;
    lim.y = ylim;
    scatter(log2(hld.lam(idx.lam{i}(idx.rho{i}))), log2(hld.rho(idx.rho{i})), 200, 'p',...
      'markeredgecolor', 'w',...
      'markerfacecolor', 'w');
    hold off;
    colormap(gca, 'parula');
    xlabel('$\lambda$',...
      'interpreter', 'latex',...
      'fontsize', 16);
    ylabel('$\rho$',...
      'interpreter', 'latex',...
      'fontsize', 16);
      set(colorbar, 'ytick', minmax(dyn.(tmp{i})));
    set(gca,...
      'xtick', tick.x,...
      'ytick', tick.y,...
      'ticklabelinterpreter', 'latex',...
      'xticklabel', ticklab.x,...
      'yticklabel', ticklab.y,...
      'fontsize', 16);
    axis([-1.05 2.05 -30.5 0.5]);
    if header.pr 
      tmp2 = sprintf('nrmse,w-%s.eps', tmp{i});
      print('-depsc', tmp2);
    end
  end
end
