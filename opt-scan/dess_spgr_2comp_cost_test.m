% script dess_spgr_2comp_cost_test.m
% header file to test dess_spgr_2comp_cost(...)
%
% copyright 2016, gopal nataraj, university of michigan
%
% version control
%   1.1     2016-04-07      original
%   1.2     2016-08-16      changed format of P

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

%% scan parameter object 
% arbitrary flip angles and repetition times
% S.de = 10;
% S.sp = 10;
% P.de.aex  = (pi/180) * linspace(1, 180, S.de)';   % rad
% P.sp.aex  = (pi/180) * linspace(1, 180, S.sp)';   % rad
% P.de.tr   = 20 * ones(S.de,1);                    % ms
% P.sp.tr   = 15 * ones(S.sp,1);                    % ms

% optimized flip angles and repetition times, from paper
P.de.aex    = (pi/180) * col([32.995 18.303 15.147]);
P.sp.aex    = [];
P.de.tr     = col([17.500 30.197 60.303]);
P.sp.tr     = [];

%% evaluate cost at different object parameter distribution sampling densities
% coarser sampling (faster)
tic;
c.test = dess_spgr_2comp_cost(P,...
  'x.ff.nsamp', 5,...
  'x.T1f.nsamp', 1,...
  'x.T1s.nsamp', 1,...
  'x.T2f.nsamp', 1,... 
  'x.T2s.nsamp', 1,...
  'x.kfs.nsamp', 1,...
  'nu.kap.nsamp', 3);
t.test = toc;
fprintf('Coarser sampling:\n');
fprintf('  Estimate of expected cost: %0.4f\n', c.test);
fprintf('  Run time: %0.3f ms\n', t.test*1000);

% finer sampling (slower, default)
tic;
c.def = dess_spgr_2comp_cost(P);
t.def = toc;
fprintf('Finer sampling:\n');
fprintf('  Estimate of expected cost: %0.4f\n', c.def);
fprintf('  Run time: %0.3f s\n', t.def);
