% script mri_mese_mwf_test.m
% test script for mri_mese_mwf(...)
% 
% gopal nataraj
% copyright 2017, university of michigan
%
% version control
%   2017-11-04      original
%   2018-02-22      adapted to test l2 regularization

% irt setup
if ~exist('irtdir', 'var')
  curdir = cd('~/Box Sync/work/irt');
  irtdir = pwd;
  setup();
  cd(curdir);
end

% add relevant directories
addpath('../../exp/sim/');
addpath('../../model/mese');
addpath('../../etc');

% constant declarations
c.b0 = 3.0;                             % field strength
c.sl = 81;                              % slice number
c.sig = 3.8607e-4;                      % noise standard deviation

c.rng.kap.tru = 0.2;
c.rng.kap.est = 0.2;

% header options
header.t1 = 0;                          % assign diff t1 to each comp but only assume bulk t1 known
header.im = 1;                          % show images
header.pr = 0;                          % print images

% 3-compartment weights
%   rows are [mw; iew; fw]
%   cols are [csf, gm, wm]
comp.wght = [0  0.03  0.15;...
             0  0.94  0.82;...
             1  0.03  0.03];
[dim.comp, dim.roi] = size(comp.wght);
           
% nominal compartmental values [mw; iew; fw]
comp.xnom.t1  = [ 500; 1000;  3000];    % ms 
comp.xnom.t2  = [  20;   80;  3000];    % ms

% load digital phantom file
f.name = 'phantom_1.0mm_msles2_crisp.fld';
if exist('flip', 'builtin')
  f.label = flip(fld_read(f.name, 'slice', c.sl), 2);
else
  f.label = flipdim(fld_read(f.name, 'slice', c.sl), 2);
end
[dim.x, dim.y] = size(f.label);

% true latent parameter maps
[xx, yy] = ndgrid(linspace(-1,1,dim.x), linspace(-1,1,dim.y));

x.tot = zeros(dim.x, dim.y);
x.mw  = zeros(dim.x, dim.y);
x.iew = zeros(dim.x, dim.y);
x.fw  = zeros(dim.x, dim.y);

nu.t1   = zeros(dim.x, dim.y);
nu.kap  = (1+c.rng.kap.tru)-c.rng.kap.tru*(xx.^2+yy.^2);

for r = 0:10
  tmp = mri_brainweb_params(r, 'b0', c.b0);
  switch r
    case {1,2,3}                        % csf, gm, wm
      x.mw(f.label==r)  = comp.wght(1,r);
      x.iew(f.label==r) = comp.wght(2,r);
      x.fw(f.label==r)  = comp.wght(3,r);
      x.tot(f.label==r) = tmp.pd;
      nu.t1(f.label==r) = tmp.t1;       % assume only bulk t1 known
  end
end

% masks
mask.wm = f.label==3;
mask.gm = f.label==2;
mask.csf = f.label==1;
mask.disp = mask.wm | mask.gm | mask.csf;
mask.est = mask.disp;
      
% scan parameters
dim.e     = 32;
P.ex.a    = pi/2;
P.ex.ph   = pi/2;
P.ref.a   = ones(dim.e,1)*pi;
P.ref.ph  = ones(dim.e,1)*0;
P.ncyc    = ones(dim.e,2);
P.tr      = 1200;
P.te      = ones(dim.e,1)*10;

% noiseless multicompartmental signals
y = zeros(dim.x, dim.y, dim.e);
tmp = fieldnames(x);
for i = 1:dim.comp
  if header.t1
    tmp1 = bsxfun(@times, comp.xnom.t1(i), ones(dim.x, dim.y));
  else
    tmp1 = nu.t1;
  end
  tmp2 = bsxfun(@times, comp.xnom.t2(i), ones(dim.x, dim.y));
  y = y + bsxfun(@times, x.(tmp{i+1}),...
    mese(x.tot, tmp1, tmp2, nu.kap, P.ex, P.ref, P.ncyc, P.tr, P.te, 'mask', mask.est));
end

% add complex gaussian noise
rng('default');
n = c.sig * (randn(size(y))+1i*randn(size(y)));
y = y + n;

% use magnitude data
y = abs(y);

% compute snr
opt.snr = {'unit', 'amp'};
snr.wm  = snr_gn(y, n, mask.wm, opt.snr{:});
snr.gm  = snr_gn(y, n, mask.gm, opt.snr{:});
snr.csf = snr_gn(y, n, mask.csf, opt.snr{:});

% parameter estimation options
dist.boxcon = [10 3000];
dist.nsamp = 100;
dist.prior = 'logunif';
reg = 2^-23;
kmean.C = 100;
bool.chat = 1;
bool.scale = 1;
bool.clust = 1;
bool.reg = 1;
bool.pool = 0;
bool.norm = 1;

% instantiate parallel workers
if bool.pool && exist('parfor', 'builtin') && isempty(gcp)
  tmp = col([2 16]);
  pool = parpool(tmp);
end

% parameter estimation
opt.map = {...
  'mask', mask,...
  'nu', nu,...
  'dist', dist,...
  'reg', reg,...
  'kmean.C', kmean.C,...
  'bool', bool
};
[xhat, t, t2] = mri_mese_mwf(y, P, opt.map{:});

if header.im
  % dynamic ranges
  dyn.tot = [0 1];
  dyn.mw  = [0 0.3];
  dyn.iew = [0 1];
  dyn.fw  = [0 1];
  
  % options
  opt.im.meth  = {'truth', 'nnls'};
  
  % compartmental maps
  tmp = fieldnames(x);
  for i = 2 %1:length(tmp)
    figure; 
    im('notick', cat(3, x.(tmp{i}), xhat.(tmp{i})), dyn.(tmp{i}), 'cbar', ' ');
    colormap(gca, 'hot');
    text(col(dim.x/2-1:dim.x:3*dim.x/2), zeros(2,1), col(opt.im.meth),...
      'VerticalAlignment', 'bottom',...
      'HorizontalAlignment', 'center',...
      'FontSize', 16,...
      'Color', 'k');
    if header.pr
      tmp = sprintf('%s,sl-%u,log2reg-n%u.eps', tmp{i}, c.sl, -log2(reg));
      print('-depsc', tmp);
    end
  end
  
  % t2 distribution map
  figure;
  im('col', dist.nsamp/10, t2.dist, 'cbar');
  colormap(gca, 'hot');
  hold on;
  count = 0;
  for i = 1:10
    for j = 1:dist.nsamp/10
      count = count+1;
      text(dim.x*j, dim.y*i, sprintf('%.1fms', t2.samp(count)),...
        'verticalalignment', 'bottom',...
        'horizontalalignment', 'right',...
        'fontsize', 6,...
        'color', 'w');
    end
  end
  hold off;
  if header.pr
    tmp = sprintf('t2-dist,sl-%u,log2reg-n%u.eps', c.sl, -log2(reg));
    print('-depsc', tmp);
  end
end
