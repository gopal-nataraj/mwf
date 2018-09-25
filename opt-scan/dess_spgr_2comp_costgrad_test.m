% script dess_spgr_2comp_costgrad_test.m
% header file to test dess_spgr_2comp_costgrad(...)
% 
% copyright 2016, gopal nataraj, university of michigan
%
% version control
%   1.1   2016-04-26      original
%   1.2   2016-04-29      exploring step-size modifications
%   1.2   2016-08-16      changed format of gradP

% setup
if ~exist('irtdir', 'var')
  curdir = cd('~/Box Sync/work/irt');
  irtdir = pwd;
  setup(); 
  cd(curdir);
end
addpath('../model/spgr/');
addpath('../model/dess/');
addpath('../crb/');
addpath('../etc/');

% construct parameter object
S.de = 10;
S.sp = 10;
TRd_rng = [17.5 30];
TRs_rng = [12.2 20];
fd_rng = [2 178] * (pi/180);
fs_rng = [2 178] * (pi/180);

P0.de.aex = linspace(fd_rng(1), fd_rng(2), S.de)';
P0.sp.aex = linspace(fs_rng(1), fs_rng(2), S.sp)';
P0.de.tr = TRd_rng(1) + (TRd_rng(2)-TRd_rng(1)) .* rand(S.de,1);
P0.sp.tr = TRs_rng(1) + (TRs_rng(2)-TRs_rng(1)) .* rand(S.sp,1);

% P0.de.aex = (pi/180) * linspace(1, 180, S.de)';   % rad
% P0.sp.aex = (pi/180) * linspace(1, 180, S.sp)';   % rad
% P0.de.tr = 20 * ones(S.de,1);                      % ms
% P0.sp.tr = 15 * ones(S.sp,1);                      % ms

% cost function parameters
costArg = {...
  'x.ff.nsamp', 5,...
  'x.T1f.nsamp', 1,...
  'x.T1s.nsamp', 1,...
  'x.T2f.nsamp', 1,...
  'x.T2s.nsamp', 1,...
  'x.kfs.nsamp', 1,...
  'nu.kap.nsamp', 3};

% gradient function parameters
gradArg = {...
  'x.ff.nsamp', 5,...
  'x.T1f.nsamp', 1,...
  'x.T1s.nsamp', 1,...
  'x.T2f.nsamp', 1,...
  'x.T2s.nsamp', 1,...
  'nu.kap.nsamp', 3};

% evaluate gradient at initial point
grad = dess_spgr_2comp_costgrad(P0, gradArg{:});

% evaluate cost after taking step in negative gradient direction
minStep = 1e-3;
maxStep = 1e1;
nStep = 50;
step = col(logspace(log10(minStep), log10(maxStep), nStep));
cost = NaN(nStep, 1);

for i = 1:nStep
  % step in negative gradient direction
  Pn.de.aex  = P0.de.aex - step(i) * grad.de.aex;
  Pn.sp.aex  = P0.sp.aex - step(i) * grad.sp.aex;
  Pn.de.tr    = P0.de.tr - step(i) * grad.de.tr;
  Pn.sp.tr    = P0.sp.tr - step(i) * grad.sp.tr;
  
  % abort early if out-of-bounds
  z.de.aex = abs(rem(Pn.de.aex, pi/2)) <= 0.01;
  z.sp.aex = abs(rem(Pn.sp.aex, pi/2)) <= 0.01;
  z.de.tr   = abs(Pn.de.tr) <= 0.01;
  z.sp.tr   = abs(Pn.sp.tr) <= 0.01;
  if sum([z.de.aex; z.sp.aex; z.de.tr; z.sp.tr]) > 0
    fprintf('gradient singularity at step = %0.3f.\n', step(i));
    break;
  end
  
  cost(i) = dess_spgr_2comp_cost(Pn, costArg{:});
end

% display
% cost descends with small steps in negative gradient direction
% cost begins to increase with larger steps due to its non-convexity
figure; semilogx(step, cost); grid on;
xlabel('step size'); 
ylabel('cost');
 