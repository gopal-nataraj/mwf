% mwf mapping in simulation with optimized dess vs multi-echo spin echo
%
% gopal nataraj
% copyright 2018, university of michigan
%
% version control
%   2018-02-24      original
%   2018-09-16      added option to only print dess-perk results for paper

% irt setup
if ~exist('irtdir', 'var')
  curdir = cd('~/Box Sync/work/irt');
  irtdir = pwd;
  setup();
  cd(curdir);
end

% add relevant directories
addpath('..');
addpath('../../../model/spgr');
addpath('../../../model/dess');
addpath('../../../model/mese');
addpath('../../../map/dess-spgr');
addpath('../../../map/mese');
addpath('../../../etc');

% constant declarations
c.b0 = 3.0;                             % field strength
c.sl = 81;                              % slice number
c.sig.dess = 3.8607e-4;                 % noise standard deviation
c.sig.mese = c.sig.dess / sqrt(2);      % mimics 2x signal averaging

c.rng.kap.tru = 0.2;
c.rng.kap.est = 0.2;
c.rng.b0f.tru  = 0.00;                  % kHz
c.rng.b0f.est  = 0.00;                  % kHz
c.rng.b0s.tru  = 0.00;                  % kHz
c.rng.b0s.est  = 0.00;                  % kHz

% header options
header.t1 = 0;                          % assign diff t1 to each comp but only assume bulk t1 known
header.im = 1;                          % show images
header.sv = 0;                          % save estimates
header.pr = 0;                          % print images to file
header.st = 1;                          % print sample statistics to file
header.ml = 1;                          % include dess-ml results

% 2-compartment weights (in gm/wm)
%   rows are [mw; iew; fw]
%   cols are [csf, gm, wm]
comp.wght = [0  0.03  0.15;...
             0  0.97  0.85;...
             1  0.00  0.00];
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
[xx, yy] = ndgrid(linspace(-1,1,dim.x), linspace(-1,1,dim.y));

% mese latent parameter maps
xm.tot = zeros(dim.x, dim.y);
xm.mw  = zeros(dim.x, dim.y);
xm.iew = zeros(dim.x, dim.y);
xm.fw  = zeros(dim.x, dim.y);

% dess-only latent parameter maps
xd.m0  = zeros(dim.x, dim.y);
xd.ff  = zeros(dim.x, dim.y);
xd.t1f = zeros(dim.x, dim.y);
xd.t1s = zeros(dim.x, dim.y);
xd.t2f = zeros(dim.x, dim.y);
xd.t2s = zeros(dim.x, dim.y);

% mese latent parameter maps
num.kap = (1+c.rng.kap.tru/2) - c.rng.kap.tru*(xx.^2 + yy.^2);
num.t1 = zeros(dim.x, dim.y);

% dess latent parameter maps
nud.kap = num.kap;
nud.b0f = (0+c.rng.b0f.tru/2) - c.rng.b0f.tru*(xx.^2 + yy.^2);
nud.b0s = (0+c.rng.b0s.tru/2) - c.rng.b0s.tru*(xx.^2 + yy.^2);
nud.r2pf = zeros(dim.x, dim.y);
nud.r2ps = zeros(dim.x, dim.y);

for r = 0:10
  tmp = mri_brainweb_params(r, 'b0', c.b0);
  switch r
    case {1,2,3}                        % csf, gm, wm
      xm.tot(f.label==r)    = tmp.pd;
      xm.mw(f.label==r)     = comp.wght(1,r);
      xm.iew(f.label==r)    = comp.wght(2,r);
      xm.fw(f.label==r)     = comp.wght(3,r);
      
      xd.m0(f.label==r)     = tmp.pd;
      xd.ff(f.label==r)     = comp.wght(1,r);
      xd.t1f(f.label == r)  = comp.xnom.t1(1);
      xd.t1s(f.label == r)  = comp.xnom.t1(2);
      xd.t2f(f.label == r)  = comp.xnom.t2(1);
      xd.t2s(f.label == r)  = comp.xnom.t2(2);
      
      num.t1(f.label==r)    = tmp.t1; 
  end
end

% masks
mask.wm = f.label==3;
mask.gm = f.label==2;
mask.csf = f.label==1;
mask.disp = mask.wm | mask.gm | mask.csf;
mask.est = mask.disp;
mask.noise = ~imdilate(mask.disp, strel('disk', 20));

% mese scan parameters
dim.e       = 32;
Pm.ex.a     = pi/2;
Pm.ex.ph    = pi/2;
Pm.ref.a    = ones(dim.e,1)*pi;
Pm.ref.ph   = ones(dim.e,1)*0;
Pm.ncyc     = ones(dim.e,2);
Pm.tr       = 600;
Pm.te       = ones(dim.e,1)*10;

% noiseless multicompartmental mese signals
sm = zeros(dim.x, dim.y, dim.e);
tmp = fieldnames(xm);
for i = 1:dim.comp
  if header.t1
    tmp1 = bsxfun(@times, comp.xnom.t1(i), ones(dim.x, dim.y));
  else
    tmp1 = num.t1;
  end
  tmp2 = bsxfun(@times, comp.xnom.t2(i), ones(dim.x, dim.y));
  sm = sm + bsxfun(@times, xm.(tmp{i+1}),...
    mese(xm.tot, tmp1, tmp2, num.kap, Pm.ex, Pm.ref, Pm.ncyc, Pm.tr, Pm.te,...
    'mask', mask.est));
end

% dess scan parameters
Pd.de.tr = [17.500; 30.197; 60.303];        % ms
Pd.de.te = ones(3,1) * [5.29 5.29];         % ms
Pd.de.aex = [32.995; 18.303; 15.147];       % deg
Pd.de.aex = Pd.de.aex * (pi/180);           % rad
Sd.de = size(Pd.de.te,1);
Ed.de = size(Pd.de.te,2);

% noiseless multicompartmental dess signals
sd = NaN(dim.x, dim.y, Sd.de, Ed.de);
for s = 1:Sd.de
  [sd(:,:,s,1), sd(:,:,s,2)] = dess_2comp(...
    xd.m0, xd.ff, xd.t1f, xd.t1s, xd.t2f, xd.t2s,...
    nud.kap, nud.b0f, nud.b0s,...
    Pd.de.aex(s), Pd.de.tr(s), Pd.de.te(s,1), Pd.de.te(s,2),...
    'R2p_f', nud.r2pf,...
    'R2p_s', nud.r2ps,...
    'exchg', 0,...
    'mag', 0);
end

% add complex gaussian noise
rng('default');
n.mese = c.sig.mese * (randn(size(sm))+1i*randn(size(sm)));
ym = sm + n.mese;
n.dess = c.sig.dess * (randn(size(sd))+1i*randn(size(sd)));
yd.de = sd + n.dess;

% use magnitude data
ym = abs(ym);
yd.de = abs(yd.de);

% compute snr
opt.snr   = {'unit', 'amp'};
opt.scan  = {'mese', 'dess'};
opt.roi   = {'wm', 'gm', 'csf'};
for i = 1:length(opt.scan)
  switch i
    case 1
      tmp = ym;
    case 2
      tmp = yd.de;
    otherwise
      error('scan unknown!?');
  end
  for j = 1:length(opt.roi)
    snr.(opt.scan{i}).(opt.roi{j}) = ...
      snr_gn(tmp, n.(opt.scan{i}), mask.(opt.roi{j}), opt.snr{:});
  end
end

% mese parameter estimation options
dist.boxcon = [10 3000];
dist.nsamp = 100;
dist.prior = 'logunif';
reg = 2^-13; 
kmean.C = 100;
bool.chat = 1;
bool.scale = 1;
bool.clust = 1;
bool.reg = [];
bool.pool = 0;
bool.norm = 1;

% instantiate parallel workers
if bool.pool && exist('parfor', 'builtin') && isempty(gcp)
  tmp = col([2 16]);
  pool = parpool(tmp);
end

% mese parameter estimation
try
  fprintf('\nTrying to load mese parameter maps...');
  load(sprintf('xhat,mese,sl-%u.mat', c.sl));
  fprintf('success!\n');
catch
  fprintf('failed: will estimate.\n');
  xhat.mese = struct('tot', [], 'mw', [], 'iew', [], 'fw', []);
  t.mese = nan(2,1);
  t2.mese = cell(2,1);
  for p = 1:2
    switch p
      case 1
        bool.reg = 0;
      case 2
        bool.reg = 1;
    end
    opt.mese = {...
      'mask', mask,...
      'nu', num,...
      'dist', dist,...
      'reg', reg,...
      'kmean.C', kmean.C,...
      'bool', bool
    };
    [xhat.mese(p), t.mese(p), t2.mese{p}] = mri_mese_mwf(ym, Pm, opt.mese{:});
  end
  
  % save xhat maps
  if header.sv
    save(sprintf('xhat,mese,sl-%u.mat', c.sl), 'xhat', 't', 't2');
  end
end

% dess parameter estimation options
meth.init = {'vpm','krr'};
meth.iter = 'pgpm';
dist.x.ff.supp = [-0.1 0.4]';
dist.x.ff.nsamp = 11;
% dist.x.ff.nsamp = 41;
dist.x.ff.prior = 'unif';
dist.x.t1f.supp = [400 600]';
dist.x.t1f.nsamp = 6;
% dist.x.t1f.nsamp = 21;
dist.x.t1f.prior = 'logunif';
dist.x.t1s.supp = [800 1200]';
dist.x.t1s.nsamp = 6;
% dist.x.t1s.nsamp = 21;
dist.x.t1s.prior = 'logunif';
dist.x.t2f.supp = [16 24]';
dist.x.t2f.nsamp = 6;
% dist.x.t2f.nsamp = 21;
dist.x.t2f.prior = 'logunif';
dist.x.t2s.supp = [64 96]';
dist.x.t2s.nsamp = 6;
% dist.x.t2s.nsamp = 21;
dist.x.t2s.prior = 'logunif';
kmean.C = 20;
rff.snr = 0.1;
rff.std = [];
rff.len = [];
rff.c = 2^0.3;                            
rff.H = 10^3;
rff.K = 10^6;
inv.krr = 2^-19;
stop.iter = [0 0];
bool.exchg  = 0;                           
bool.mag.sp = 1;
bool.mag.de = 1;
bool.chat = 1;
bool.norm = 1; 
bool.reset = 1;
bool.rfftst = 0;
bool.nuclip = 1;
bool.reg = 0;
bool.precon = 1;
bool.disp = 0;

% dess parameter estimation
try
  fprintf('\nTrying to load dess parameter maps...');
  tmp = xhat;
  tmp2 = t;
  load(sprintf('xhat,dess,sl-%u.mat', c.sl));
  xhat = catstruct(tmp, xhat);
  t = catstruct(tmp2, t);
  fprintf('success!\n');
catch
  fprintf('failed: will estimate.\n');
  xhat.dess = struct('init', [], 'iter', []);
  t.dess = struct('init', [], 'iter', []);
  for p = 1:length(meth.init)
    if strcmp(meth.init(p), 'vpm')
      opt.map = {...
        'mask.disp', mask.disp,...
        'mask.est', mask.est,...
        'nu', nud,...
        'meth.init', meth.init{p},...
        'meth.iter', meth.iter,...
        'dist.x.ff', dist.x.ff,...
        'dist.x.t1f', dist.x.t1f,...
        'dist.x.t1s', dist.x.t1s,...
        'dist.x.t2f', dist.x.t2f,...
        'dist.x.t2s', dist.x.t2s,...
        'kmean.C', kmean.C,...
        'stop.iter', stop.iter(p),...
        'bool', bool...
      };
    elseif strcmp(meth.init(p), 'krr')
      opt.map = {...
        'mask.disp', mask.disp,...
        'mask.est', mask.est,...
        'mask.noise', mask.noise,...
        'nu', nud,...
        'meth.init', meth.init{p},...
        'kmean.C', kmean.C,...
        'rff', rff,...
        'inv.krr', inv.krr,...
        'stop.iter', stop.iter(p),...
        'bool', bool...
      };
    else
      error('unknown meth.init!');
    end
    [xhat.dess(p), t.dess(p)] = mri_multicomp_map(yd, Pd, opt.map{:});
  end
  
  % save xhat maps
  if header.sv
    tmp = xhat;
    xhat = rmfield(xhat, {'mese'});
    tmp2 = t;
    t = rmfield(t, {'mese'});
    save(sprintf('xhat,dess,sl-%u.mat', c.sl), 'xhat', 't');
    xhat = tmp;
    t = tmp2;
  end
end

if header.im
  % mese dynamic ranges
  dyn.mese.tot  = [0 1];
  dyn.mese.mw   = [0 0.3];
  dyn.mese.iew  = [0 1];
  dyn.mese.fw   = [0 1];
  
  % dess dynamic ranges
  dyn.dess.m0   = [0 1];
  dyn.dess.ff   = [0 0.3];
  dyn.dess.t1f  = [0 1000];
  dyn.dess.t1s  = [0 2000];
  dyn.dess.t2f  = [0 50];
  dyn.dess.t2s  = [0 200];
  
  % options
  opt.im.mese = {'Truth','NNLS','RNNLS'};
  opt.im.dess = {'Truth','ML','PERK'};
  if header.ml
    opt.im.both = {'Truth','MESE-NNLS','MESE-RNNLS','DESS-ML','DESS-PERK'};
  else
    opt.im.both = {'Truth','MESE-NNLS','MESE-RNNLS','DESS-PERK'};
  end
  opt.unit.mese = {'a.u.',' ',' ',' '};
  opt.unit.dess = {'a.u.',' ','ms','ms','ms','ms'};
  opt.type = {'im','err'};
  opt.roi = {'wm','gm'};
  
  % mese-only maps
  tmp = fieldnames(xhat.mese);
  for i = 1:length(tmp)
    for j = 1:length(opt.type)
      tmp2 = cat(3, xm.(tmp{i}), xhat.mese(1).(tmp{i}), xhat.mese(2).(tmp{i}));
      if strcmp(opt.type{j}, 'err')
        tmp2 = bsxfun(@minus, tmp2, xm.(tmp{i}));
        tmp2 = abs(tmp2);
      end
      tmp2 = embed(masker(tmp2, mask.wm|mask.gm), mask.wm|mask.gm);
      figure; im('notick', tmp2, dyn.mese.(tmp{i}), 'cbar', ' ');
      colormap(gca, 'hot');
      tmp2 = colorbar;
      ylabel(tmp2, opt.unit.mese{i});
      set(tmp2, 'ytick',...
        linspace(dyn.mese.(tmp{i})(1), dyn.mese.(tmp{i})(2), 6));
      hold on;
      if strcmp(opt.type{j},'im')
        text(dim.x/2-1:dim.x:5*dim.x/2, zeros(3,1), col(opt.im.mese),...
          'VerticalAlignment', 'bottom',...
          'HorizontalAlignment', 'center',...
          'FontSize', 16,...
          'Color', 'k');
      end
      if strcmp(opt.type{j}, 'im')
        tmp2 = upper(tmp{i});
      else
        tmp2 = [upper(tmp{i}) ' Magnitude Error'];
      end
      text(0, dim.y/2-1, tmp2,...
        'VerticalAlignment', 'bottom',...
        'HorizontalAlignment', 'center',...
        'FontSize', 14,...
        'Color', 'k',...
        'Rotation', 90);
      hold off;
      if header.pr
        tmp2 = sprintf('mese-%s,sl-%u,%s.eps',...
          tmp{i}, c.sl, opt.type{j});
        print('-depsc', tmp2);
      end
    end
  end
  
  % dess-only maps
  tmp = fieldnames(xhat.dess(1).init);
  for i = 1:length(tmp)
    for j = 1:length(opt.type)
      tmp2 = cat(3, xd.(tmp{i}), xhat.dess(1).init.(tmp{i}), xhat.dess(2).init.(tmp{i}));
      if strcmp(opt.type{j}, 'err')
        tmp2 = bsxfun(@minus, tmp2, xd.(tmp{i}));
        tmp2 = abs(tmp2);
      end
      tmp2 = embed(masker(tmp2, mask.wm|mask.gm), mask.wm|mask.gm);
      figure; im('notick', tmp2, dyn.dess.(tmp{i}), 'cbar', ' ');
      colormap(gca, 'hot');
      tmp2 = colorbar;
      ylabel(tmp2, opt.unit.dess{i});
      set(tmp2, 'ytick',...
        linspace(dyn.dess.(tmp{i})(1), dyn.dess.(tmp{i})(2), 6));
      hold on;
      if strcmp(opt.type{j},'im')
        text(dim.x/2-1:dim.x:5*dim.x/2, zeros(3,1), col(opt.im.dess),...
          'VerticalAlignment', 'bottom',...
          'HorizontalAlignment', 'center',...
          'FontSize', 16,...
          'Color', 'k');
      end
      if strcmp(opt.type{j}, 'im')
        tmp2 = upper(tmp{i});
      else
        tmp2 = [upper(tmp{i}) ' Magnitude Error'];
      end
      text(0, dim.y/2-1, tmp2,...
        'VerticalAlignment', 'bottom',...
        'HorizontalAlignment', 'center',...
        'FontSize', 14,...
        'Color', 'k',...
        'Rotation', 90);
      hold off;
      if header.pr
        tmp2 = sprintf('dess-%s,sl-%u,%s.eps',...
          tmp{i}, c.sl, opt.type{j});
        print('-depsc', tmp2);
      end
    end
  end
  
  % mese-mwf, dess-ff comparison
  for j = 1:length(opt.type)
    if header.ml
      tmp = cat(3, xm.mw,...
        xhat.mese(1).mw, xhat.mese(2).mw,...
        xhat.dess(1).init.ff, xhat.dess(2).init.ff);
    else
      tmp = cat(3, xm.mw,...
        xhat.mese(1).mw, xhat.mese(2).mw,...
        xhat.dess(2).init.ff);
    end
    if strcmp(opt.type{j}, 'err')
      tmp = bsxfun(@minus, tmp, xm.mw);
      tmp = abs(tmp);
    end
    tmp = embed(masker(tmp, mask.wm|mask.gm), mask.wm|mask.gm);
    figure; im('notick', 'row', 1, tmp, dyn.mese.mw, 'cbar', ' ');
    colormap(gca, 'hot');
    tmp = colorbar;
    set(tmp, 'ytick',...
      linspace(dyn.mese.mw(1), dyn.mese.mw(2), 6));
    hold on;
    if strcmp(opt.type{j}, 'im')
      tmp = length(opt.im.both);
      text(dim.x/2-1:dim.x:(tmp-1/2)*dim.x, zeros(tmp,1), col(opt.im.both),...
        'VerticalAlignment', 'bottom',...
        'HorizontalAlignment', 'center',...
        'FontSize', 11,...
        'Color', 'k');
    end
    if strcmp(opt.type{j}, 'im')
      tmp = 'MWF';
    else
      tmp = 'MWF Mag. Error';
    end
    text(0, dim.y/2-1, tmp,...
      'VerticalAlignment', 'bottom',...
      'HorizontalAlignment', 'center',...
      'FontSize', 12,...
      'Color', 'k',...
      'Rotation', 90);
    hold off;
    if header.pr
      tmp = sprintf('mese-mw,dess-ff,sl-%u,%s.eps', c.sl, opt.type{j});
      print('-depsc', tmp);
    end
  end
end

% mese-mwf, dess-ff summary statistics
if header.st
  tmp = sprintf('mese-mw,dess-ff,sl-%u,stat', c.sl);
  fid = fopen(tmp, 'w');
  fprintf(fid, 'mw/ff estimator statistics for slice %u\n', c.sl);
  fprintf(fid, '\trows denote regions of interest\n');
  fprintf(fid, '\tcols denote acquisitions-estimators\n');  

  tmp = '';
  for i = 1:length(opt.im.both)
    tmp = strcat(tmp, '%28s');
  end
  tmp = strcat(tmp, '\n');
  fprintf(fid, tmp, opt.im.both{:});
end
for r = 1:length(opt.roi)
  for m = 1:length(opt.im.both)
    switch opt.im.both{m}
      case 'Truth'
        tmp = xm.mw;
      case 'MESE-NNLS'
        tmp = xhat.mese(1).mw;
      case 'MESE-RNNLS'
        tmp = xhat.mese(2).mw;
      case 'DESS-ML'
        tmp = xhat.dess(1).init.ff;
      case 'DESS-PERK'
        tmp = xhat.dess(2).init.ff;
    end
    summ(m,r) = stat(...
      masker(tmp, mask.(opt.roi{r})),...
      'true', mean(masker(xm.mw, mask.(opt.roi{r}))));
    if header.st
      fprintf(fid, '\t%7.4f\t%c%7.4f (%6.4f)',...
        summ(m,r).mean,...
        char(177),...
        summ(m,r).std,...
        summ(m,r).rmse);
    end
  end
  if header.st
    fprintf(fid, '\n');
  end
end
if header.st
  fclose(fid);
end
