% mwf mapping in vivo with optimized dess vs multi-echo spin echo
%
% written by: gopal nataraj
% copyright 2017, university of michigan
%
% version
%   2018-01-24    acquired data
%   2018-01-25    original script, copied from 2018-01-10
%   2018-02-26    added roi selection, pretty plots

% irt
if ~exist('irtdir', 'var')
  curdir = cd('~/Box Sync/work/irt');
  irtdir = pwd;
  setup();
  cd(curdir);
end

% mapping
addpath('../../../../../relax/mat/map/b1');
addpath('../../../../../relax/mat/map/sense');
addpath('../../../../../relax/mat/map/t1-t2');
addpath('../../../../../relax/mat/model/spgr');
addpath('../../../map');
addpath('../../../model/mese');
addpath('../../../model/dess');
addpath('../../../etc');

% header options
header.sv   = 0;                          % save coil-comb, param reconstructions
header.rg   = 1;                          % register coil-comb images
header.rs   = 0;                          % roi manual selection
header.im   = 1;                          % show images
header.roi  = 1;                          % show rois in images
header.pr   = 0;                          % print images to file
header.stat = 1;                          % print sample statistics to file

% constant declarations
dim.o.lo = [200 50 8];               
dim.o.hi = [200 200 8];

dim.c = 32;
dim.sl = 5;
dim.e = 32;
dim.r = 2;                                % number of mese repetitions

dim.fov = [240 240 24];                   % mm

% load image data
load(sprintf('im,kap,sl-%u.mat', dim.sl), 'nu');
tmp = nu;
load(sprintf('im,t1,sl-%u.mat', dim.sl), 'nu');
nu = catstruct(tmp, nu);
load(sprintf('im,mese,sl-%u.mat', dim.sl), 'ysos', 'ycc');
tmp = ysos;
tmp2 = ycc;
load(sprintf('im,dess,sl-%u.mat', dim.sl), 'ysos', 'ycc');
ysos = catstruct(tmp, ysos);
ycc = catstruct(tmp2, ycc);

% % register coil-combined images and nu.* to first mese image
if header.rg 
  fprintf('\nRegistering coil-combined images to each other:\n');
  ysos.mese = permute(ysos.mese, [1 2 4 3]);
  ycc.mese = permute(ycc.mese, [1 2 4 3]); 
  reg.fixed = abs(ycc.mese(:,:,1,1));
  reg.ttype = 'rigid';
  [reg.optim, reg.metric] = imregconfig('multimodal'); 
  
  tmp = fieldnames(ycc);
  for i = 1:length(tmp)
    for j = 1:size(ycc.(tmp{i}),3)
      fprintf('  registering %s dataset %u...', upper(tmp{i}), j)
      tmp2 = imregtform(abs(ycc.(tmp{i})(:,:,j,1)),...
        reg.fixed, reg.ttype, reg.optim, reg.metric);
      for k = 1:size(ycc.(tmp{i}),4)
        ysos.(tmp{i})(:,:,j,k) = imwarp(ysos.(tmp{i})(:,:,j,k), tmp2,...
          'OutputView', imref2d(size(reg.fixed)));
        ycc.(tmp{i})(:,:,j,k) = imwarp(ycc.(tmp{i})(:,:,j,k), tmp2,...
          'OutputView', imref2d(size(reg.fixed)));
      end
      fprintf('done.\n');
    end
  end
  ysos.mese = permute(ysos.mese, [1 2 4 3]);
  ycc.mese = permute(ycc.mese, [1 2 4 3]); 
  
  fprintf('\nRegistering known parameter maps to coil-combined images:\n');
  tmp = fieldnames(nu);
  for i = 1:length(tmp)
    tmp2 = fieldnames(nu.(tmp{i}));
    tmp3 = imregtform(nu.(tmp{i}).(tmp2{1}), reg.fixed, reg.ttype, reg.optim, reg.metric);
    fprintf('  registering nu.%s...', tmp{i});
    for j = 1:length(tmp2)
      nu.(tmp{i}).(tmp2{1}) = imwarp(nu.(tmp{i}).(tmp2{1}), tmp3,...
        'OutputView', imref2d(size(reg.fixed)));
    end
    fprintf('done.\n');
  end 
end

% average mese repetitions
ysos.mese = mean(ysos.mese, 4);
ycc.mese = mean(ycc.mese, 4);

% create masks from mese first echo average image
tmp = squeeze(abs(ycc.mese(:,:,1)));
mask.thresh = 0.07;
mask.t = imfill(tmp > mask.thresh * max(col(tmp)), 'holes');
mask.t = imdilate(~imdilate(~mask.t, strel('disk', 5)), strel('disk', 5));
mask.b = imdilate(mask.t, strel('disk', 0));
mask.n = ~imdilate(mask.t, strel('disk', 35));

%% mese mwf estimation
% mese scan parameters
P.mese.ex.a     = pi/2;
P.mese.ex.ph    = 0;
P.mese.ref.a    = ones(dim.e,1)*pi;
P.mese.ref.ph   = ones(dim.e,1)*(pi/2);
P.mese.ncyc     = ones(dim.e,1)*[2 2];
P.mese.tr       = 600;
P.mese.te       = ones(dim.e,1)*9.20;         

% mese mwf estimation options
wght = ones(dim.e,1);
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

% mese mwf estimation
try
  fprintf('\nTrying to load mese parameter maps...');
  load(sprintf('xhat,mese,sl-%u.mat', dim.sl));
  fprintf('success!\n');
catch
  fprintf('failed: will estimate.\n');
  xhat.mese = struct('tot', [], 'mw', [], 'iew', [], 'fw', []);
  t.mese = nan(2,1);
  t2.mese = cell(2,1);
  
  % instantiate parallel workers
  if bool.pool && exist('parfor', 'builtin') && isempty(gcp)
    tmp = col([2 16]);
    pool = parpool(tmp);
  end
  
  % run without, then with l2 regularization
  for p = 1:2
    % set options
    switch p
      case 1
        bool.reg = 0;
      case 2
        bool.reg = 1;
    end
    opt.mese = {...
      'mask.disp', mask.t,...
      'mask.est', mask.b,...
      'nu.t1', nu.t1.pgpm,...
      'nu.kap', double(nu.kap.rls),...
      'wght', wght,...
      'dist', dist,...
      'reg', reg,...
      'kmean.C', kmean.C,...
      'bool', bool
    };

    % map, using magnitude echo images 
    [xhat.mese(p), t.mese(p), t2.mese{p}] = ...
      mri_mese_mwf(abs(ycc.mese), P.mese, opt.mese{:});
  end
  
  % save maps
  if header.sv
    save(sprintf('xhat,mese,sl-%u.mat', dim.sl), 'xhat', 't', 't2');
  end
end

%% dess ff estimation
% dess scan parameters
P.de.tr = [17.500; 30.197; 60.303];       % ms
P.de.te = ones(3,1) * [5.29 5.29];        % ms
P.de.aex = [32.995; 18.303; 15.147];      % deg
P.de.aex = P.de.aex * (pi/180);           % rad

% dess ff estimation options
meth.init = 'krr';
meth.iter = 'pgpm';
kmean.C = 50;
rff.snr = 0.1;
rff.std = [];
rff.len = []; 
rff.c = 2^0.3;
rff.H = 10^3;
rff.K = 10^6;
inv.krr = 2^-19;
stop.iter = 0;
stop.tolx = 10^-8;
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

% use magnitude data as appropriate
if bool.mag.de
  ysos.de = abs(ysos.de);
  ycc.de = abs(ycc.de);
end

% dess ff estimation
try
  fprintf('\nTrying to load dess parameter maps...');
  tmp = xhat;
  tmp2 = t;
  load(sprintf('xhat,dess,sl-%u.mat', dim.sl));
  xhat = catstruct(tmp, xhat);
  t = catstruct(tmp2, t);
  fprintf('success!\n');
catch
  fprintf('failed: will estimate.\n');
  
  % set options
  opt.dess = {...
    'mask.disp', mask.t,...
    'mask.est', mask.b,...
    'mask.noise', mask.n,...
    'nu.kap', double(nu.kap.rls),...
    'meth', meth,...
    'kmean.C', kmean.C,...
    'rff', rff,...
    'inv.krr', inv.krr,...
    'stop.iter', stop.iter,...
    'stop.tolx', stop.tolx,...
    'bool', bool...
  };

  % map, using magnitude sum-of-squares images for proper setting of rff.std
  tmp4.de = ysos.de;
  tmp5.de = P.de;
  [xhat.dess, t.dess] = ...
    mri_multicomp_map(tmp4, tmp5, opt.dess{:});

  % save maps
  if header.sv
    tmp = xhat;
    tmp2 = t;
    xhat = rmfield(xhat, {'mese'});
    t = rmfield(t, {'mese'});
    save(sprintf('xhat,dess,sl-%u.mat', dim.sl), 'xhat', 't');
    xhat = tmp;
    t = tmp2;
  end
end

% roi selection
if header.rs
  tmp = 'ROI selection: [m]anual or [s]aved? ';
  roi.sel = input(tmp, 's');
else
  roi.sel = 's';
end
roi.label = {'wm-ra'; 'wm-la'; 'wm-rp'; 'wm-lp';'wm-cap';'gm-cort'};
dim.roi = length(roi.label);
switch roi.sel
  case 'm'
    roi.mask = false([dim.o.hi(1:2) dim.roi]);
    for r = 1:dim.roi
      switch roi.label{r}
        case {'wm-ra','wm-la','wm-rp','wm-lp'}
          tmp = transpose(xhat.dess.iter.ff);
        case {'wm-cap'}
          tmp = transpose(xhat.mese(2).mw);
        case {'gm-cort'}
          tmp = abs(ycc.mese(:,:,1));
          tmp = bsxfun(@min, tmp, 6);
          tmp = bsxfun(@rdivide, tmp, max(col(tmp)));
          tmp = transpose(tmp);
        otherwise
          error('Unknown roi!?');
      end
      while true
        fprintf('Select %s polygon.\n', roi.label{r});
        tmp2 = roipoly(tmp);
        roi.mask(:,:,r) = roi.mask(:,:,r) | transpose(tmp2);
        tmp2 = input('Add another polygon to this ROI [*/n]? ', 's');
        if strcmp(tmp2,'n'), break; end
      end
    end
    if exist('roi-mask.mat', 'file')
      tmp = input('Overwrite previously-defined ROI mask [y/n]? ', 's');
      switch tmp
        case 'y'
          save('roi-mask.mat', 'roi'); 
          fprintf('roi-mask.mat overwritten.\n');
        case 'n'
          fprintf('New ROI mask not saved.\n');
        otherwise
          error('Unrecognized input.');
      end
    else
      save('roi-mask.mat', 'roi');
    end
  case 's'
    try
      load('roi-mask.mat');
    catch
      error('ROI masks not found!');
    end
  otherwise
    error('Unrecognized input!');
end

% roi boundaries
roi.bound = cell(dim.roi,1);
for r = 1:dim.roi
  roi.bound{r} = bwboundaries(roi.mask(:,:,r), 'noholes');
end

if header.im
  % mese dynamic ranges
  dyn.mese.tot  = [0 dim.c];
  dyn.mese.mw   = [0 0.3];
  dyn.mese.iew  = [0 1];
  dyn.mese.fw   = [0 1];
  
  % dess dynamic ranges
  dyn.dess.m0   = [0 dim.c];
  dyn.dess.ff   = [0 0.3];
  dyn.dess.t1f  = [0 1000];
  dyn.dess.t1s  = [0 2000];
  dyn.dess.t2f  = [0 50];
  dyn.dess.t2s  = [0 200];
  
  % options
  opt.scan      = {'mese','dess'};
  opt.mese      = {'NNLS','RNNLS'};
  opt.dess      = {'PERK','PERK-PGPM'};
  opt.both      = {'MESE-NNLS','MESE-RNNLS','DESS-PERK'};
  opt.x.mese    = fieldnames(dyn.mese);
  opt.x.dess    = fieldnames(dyn.dess);
  opt.unit.mese = {'a.u.',' ',' ',' '};
  opt.unit.dess = {'a.u.',' ','ms','ms','ms','ms'};
  opt.color.roi = {...
    [0.6350, 0.0780, 0.1840],...
    [0.8500, 0.3250, 0.0980],...
    [0.9290, 0.6940, 0.1250],...
    [0.4940, 0.1840, 0.5560],...
    [0.4660, 0.6740, 0.1880],...
    [0.3010, 0.7450, 0.9330],...
    [0, 0.4470, 0.7410]...
  };
  opt.color.est = {...
    [0, 0.4470, 0.7410],...
    [0.8500, 0.3250, 0.0980]...
  };
  
  % crop for display
  crop.x = 25:180;
  crop.y = 1:196;
  dim.xc = length(crop.x);
  dim.yc = length(crop.y);
    
  % mese-only and dess-only maps
  for i = 1:length(opt.scan)
    for j = 1:length(opt.x.(opt.scan{i}))
      switch opt.scan{i}
        case 'mese'
          tmp = cat(3,...
            xhat.mese(1).(opt.x.mese{j})(crop.x,crop.y),...
            xhat.mese(2).(opt.x.mese{j})(crop.x,crop.y));
        case 'dess'
          tmp = cat(3,...
            xhat.dess.init.(opt.x.dess{j})(crop.x,crop.y),...
            xhat.dess.iter.(opt.x.dess{j})(crop.x,crop.y));
      end
      figure;
      tmp2 = dyn.(opt.scan{i}).(opt.x.(opt.scan{i}){j});
      im('notick', 'row', 1, tmp, tmp2, 'cbar', ' ');
      if j==1
        colormap(gca, 'gray');
      else
        colormap(gca, 'hot');
      end
      tmp = colorbar;
      ylabel(tmp, opt.unit.(opt.scan{i}){j});
      set(tmp, 'ytick', linspace(tmp2(1), tmp2(2), 6));
      hold on;
      tmp = length(opt.(opt.scan{i}));
      text(col(dim.xc/2-1:dim.xc:(2*tmp-1)*dim.xc/2), zeros(tmp,1), col(opt.(opt.scan{i})),...
        'VerticalAlignment', 'bottom',...
        'HorizontalAlignment', 'center',...
        'FontSize', 16,...
        'Color', 'k');
      text(0, dim.yc/2-1, upper(opt.x.(opt.scan{i}){j}),...
        'VerticalAlignment', 'bottom',...
        'HorizontalAlignment', 'center',...
        'FontSize', 16,...
        'Color', 'k',...
        'Rotation', 90);
      if header.roi && j==1
        for r = 1:dim.roi
          tmp = roi.bound{r};
          for k = 1:length(roi.bound{r})
            tmp2 = tmp{k};
            tmp2 = bsxfun(@minus, tmp2, [crop.x(1)-1, crop.y(1)-1]);
            plot(tmp2(:,1), tmp2(:,2),...
              'Color', opt.color.roi{r},...
              'LineWidth', 0.5);
          end
        end
        text(dim.xc/2, 0, 'A',...
          'VerticalAlignment', 'top',...
          'HorizontalAlignment', 'center',...
          'FontSize', 16,...
          'Color', 'w');
        text(dim.xc, dim.yc/2, 'L',...
          'VerticalAlignment', 'middle',...
          'HorizontalAlignment', 'right',...
          'FontSize', 16,...
          'Color', 'w');
        text(dim.xc/2, dim.yc, 'P',...
          'VerticalAlignment', 'bottom',...
          'HorizontalAlignment', 'center',...
          'FontSize', 16,...
          'Color', 'w');
        text(1, dim.yc/2, 'R',...
          'VerticalAlignment', 'middle',...
          'HorizontalAlignment', 'left',...
          'FontSize', 16,...
          'Color', 'w');
      end
      hold off;
      if header.pr
        tmp = sprintf('%s-%s,sl-%u.eps', opt.scan{i}, opt.x.(opt.scan{i}){j}, dim.sl);
        print('-depsc', tmp);
      end
    end
  end
  
  % mese-mwf, dess-ff comparison
  tmp = cat(3,...
    xhat.mese(1).mw(crop.x,crop.y),...
    xhat.mese(2).mw(crop.x,crop.y),...
    xhat.dess.init.ff(crop.x,crop.y));
  figure;
  im('notick', 'row', 1, tmp, dyn.mese.mw, 'cbar', ' ');
  colormap(gca, 'hot');
  tmp = colorbar;
  ylabel(tmp, opt.unit.mese{2});
  set(tmp, 'ytick', linspace(dyn.mese.mw(1), dyn.mese.mw(2), 6));
  hold on;
  tmp = length(opt.both);
  text(col(dim.xc/2-1:dim.xc:(2*tmp-1)*dim.xc/2), zeros(tmp,1), col(opt.both),...
    'VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'center',...
    'FontSize', 16,...
    'Color', 'k');
  text(0, dim.yc/2-1, 'MWF',...
    'VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'center',...
    'FontSize', 16,...
    'Color', 'k',...
    'Rotation', 90);
  hold off;
  if header.pr
    tmp = sprintf('mese-mw,dess-ff,sl-%u.eps', dim.sl);
    print('-depsc', tmp);
  end

  % roi labels on mese first echo image
  if header.roi
    figure;
    im('notick', ysos.mese(crop.x,crop.y,1), [0 6], ' ');
    colormap(gca, 'gray');
    hold on;
    for r = 1:dim.roi
      tmp = roi.bound{r};
      for k = 1:length(roi.bound{r})
        tmp2 = tmp{k};
        tmp2 = bsxfun(@minus, tmp2, [crop.x(1)-1, crop.y(1)-1]);
        plot(tmp2(:,1), tmp2(:,2),...
          'Color', opt.color.roi{r},...
          'LineWidth', 0.5);
      end
    end
    text(dim.xc/2, 0, 'A',...
      'VerticalAlignment', 'top',...
      'HorizontalAlignment', 'center',...
      'FontSize', 16,...
      'Color', 'w');
    text(dim.xc, dim.yc/2, 'L',...
      'VerticalAlignment', 'middle',...
      'HorizontalAlignment', 'right',...
      'FontSize', 16,...
      'Color', 'w');
    text(dim.xc/2, dim.yc, 'P',...
      'VerticalAlignment', 'bottom',...
      'HorizontalAlignment', 'center',...
      'FontSize', 16,...
      'Color', 'w');
    text(1, dim.yc/2, 'R',...
      'VerticalAlignment', 'middle',...
      'HorizontalAlignment', 'left',...
      'FontSize', 16,...
      'Color', 'w');
    hold off;
    if header.pr
      print('-depsc', 'roi.eps');
    end
  end
end

% add a lateral-wm roi
tmp = sum(roi.mask(:,:,1:4),3)>0;
roi.mask = cat(3, tmp, roi.mask);
roi.label = {'wm-lat', roi.label{:}};
dim.roi = size(roi.mask,3);

% plot mese t2 distributions averaged over rois
for j = 1:dim.roi
  figure;
  hold on;
  for i = 1:2
    tmp = masker(t2.mese{i}.dist, stackpick(roi.mask,j));
    plot(t2.mese{i}.samp, col(mean(tmp,1)),...
      'Color', opt.color.est{i});
  end
  hold off;
  set(gca, 'XScale', 'log');
  title(roi.label{j});
  xlabel('t2 (ms)');
  legend({'NNLS', 'RNNLS'});
  if header.pr
    tmp = ['mese-dist,', roi.label{j}, sprintf(',sl-%u.eps', dim.sl)];
    print('-depsc', tmp);
  end
end

% mese-mwf, dess-ff summary statistics
if header.stat
  tmp = sprintf('mese-mw,dess-ff,sl-%u,stat', dim.sl);
  fid = fopen(tmp,'w');
  fprintf(fid, 'mw/ff estimator statistics for slice %u\n', dim.sl);
  fprintf(fid, '\trows denote regions of interest\n');
  fprintf(fid, '\tcols denote acquisitions-estimators\n');  

  tmp = '        ';
  for i = 1:length(opt.both)
    tmp = strcat(tmp, '%20s');
  end
  tmp = strcat(tmp, '\n');
  fprintf(fid, tmp, opt.both{:});
end
for r = 2:dim.roi % omit pooled wm-lat roi for consistency with paper
  if header.stat
    fprintf(fid, '%8s', roi.label{r});
  end
  for m = 1:length(opt.both)
    switch opt.both{m}
      case 'MESE-NNLS'
        tmp = xhat.mese(1).mw;
      case 'MESE-RNNLS'
        tmp = xhat.mese(2).mw;
      case 'DESS-PERK'
        tmp = xhat.dess.init.ff;
    end
    summ(m,r) = stat(masker(tmp, stackpick(roi.mask,r)));
    if header.stat
      fprintf(fid, '\t%7.4f\t%c%7.4f',...
        summ(m,r).mean,...
        char(177),...
        summ(m,r).std);
    end
  end
  if header.stat
    fprintf(fid, '\n');
  end
end
if header.stat
  fclose(fid);
end
