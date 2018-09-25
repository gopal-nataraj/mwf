% script dess_spgr_2comp_Popt_test.m
% header file to test dess_spgr_2comp_Popt(...)
%
% copyright 2016, gopal nataraj, university of michigan
%
% version control
%   1.1     2016-04-08      original
%   1.2     2016-04-20      added more options to test fmincon
%   2.1     2016-05-02      retesting when gradient also provided
%   2.2     2016-08-16      changed format of P 

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

% construct parameter initialization
C.de = 10;
C.sp = 10;
TRd_rng = [17.50 50];
TRs_rng = [12.25 30];
fd_rng = [5 90] * (pi/180);
fs_rng = [5 90] * (pi/180);

P0.de.aex = linspace(fd_rng(1), fd_rng(2), C.de)';
P0.sp.aex = linspace(fs_rng(1), fs_rng(2), C.sp)';
% P0.de.tr = TRd_rng(1) + (TRd_rng(2)-TRd_rng(1)) .* rand(C.de,1);
% P0.sp.tr = TRs_rng(1) + (TRs_rng(2)-TRs_rng(1)) .* rand(C.sp,1);
P0.de.tr = TRd_rng(1) * ones(C.de,1);
P0.sp.tr = TRs_rng(1) * ones(C.sp,1);
% P0.de.tr = 20 * ones(C.de,1);
% P0.sp.tr = 15 * ones(C.sp,1);

% linear constraints
lincon.tr = C.de*20 + C.sp*15;

% cost function options
costOpt = {...
  'x.ff.minmax', [0.03 0.21],...
  'x.ff.nsamp', 5,...
  'x.T1f.nsamp', 1,...
  'x.T1s.nsamp', 1,...
  'x.T2f.nsamp', 1,...
  'x.T2s.nsamp', 1,...
  'x.kfs.nsamp', 1,...
  'nu.kap.nsamp', 3};
  
% gradient function options
gradOpt = {...
  'x.ff.minmax', [0.03 0.21],...
  'x.ff.nsamp', 5,...
  'x.T1f.nsamp', 1,...
  'x.T1s.nsamp', 1,...
  'x.T2f.nsamp', 1,...
  'x.T2s.nsamp', 1,...
  'nu.kap.nsamp', 3};
  
% fmincon options
fminconOpt = {...
  'boxcon.de.tr', col(TRd_rng),...
  'boxcon.sp.tr', col(TRs_rng),...
  'boxcon.de.aex', col(fd_rng),...
  'boxcon.sp.aex', col(fs_rng),...
  'lincon.tr', lincon.tr,...
  'fmincon.wgrad', 'on',...
  'fmincon.alg', 'active-set',...
  'fmincon.tolFun', 1e-10,...
  'fmincon.tolX', 1e-10,...
  'fmincon.maxIter', 1000};

% parameter optimization
tic;
[P, fopt, flag] = dess_spgr_2comp_Popt(P0, costOpt, gradOpt, fminconOpt{:});
t = toc;

% compare initial and final coefficient of variation w.r.t. mean(ff)
f0        = dess_spgr_2comp_cost(P0, costOpt{:});
cfvar0    = sqrt(f0)    ./ mean([0.03 0.21]);
cfvaropt  = sqrt(fopt)  ./ mean([0.03 0.21]);
