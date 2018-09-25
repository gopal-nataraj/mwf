% script dess_2comp_gradx_hessx_test.m
% check two-compartment dess signal gradients and hessians
%
% copyright 2018, gopal nataraj, university of michigan
%
% version control
%   2018-01-16      original

% irt setup
if ~exist('irtdir', 'var')
  curdir = cd('~/Box Sync/work/irt');
  irtdir = pwd;
  setup(); 
  cd(curdir);
end

% constants
n = 1000;
bool.exchg = 0;
bool.mag = 1;

% default latent parameter values
def.x.M0 = 1;               
def.x.ff = 0.15;             
def.x.T1f = 400;              % ms
def.x.T1s = 1000;             % ms
def.x.T2f = 20;               % ms
def.x.T2s = 80;               % ms

% non-default latent parameter ranges
rng.x.M0 = [eps 1]';
rng.x.ff = [eps 1]';
rng.x.T1f = [eps 500]';       % ms
rng.x.T1s = [eps 500]';       % ms
rng.x.T2f = [eps 100]';       % ms
rng.x.T2s = [eps 100]';       % ms

% known parameter values
def.nu.kap = ones(n,1) * 1;
def.nu.Dwf = ones(n,1) * 0;   % kHz
def.nu.Dws = ones(n,1) * 0;   % kHz

% scan parameters
flip = 30 * (pi/180);         % rad
TR = 20;                      % ms
TEp = 5;                      % ms
TEm = 5;                      % ms

% loop over latent parameters
xf = fieldnames(def.x);
unit = {'a.u.',' ','ms','ms','ms','ms'};
clear('tmp');
for i = 1:length(xf)
  % set latent parameter ranges, varying one only
  for j = 1:length(xf)
    if i==j
      tmp.(xf{j}) = col(linspace(rng.x.(xf{j})(1), rng.x.(xf{j})(2), n));
    else
      tmp.(xf{j}) = ones(n,1) * def.x.(xf{j});
    end
  end
  
  % 2-compartment dess signals
  [sp, sm] = dess_2comp(...
    tmp.M0, tmp.ff, tmp.T1f, tmp.T1s, tmp.T2f, tmp.T2s,...
    def.nu.kap, def.nu.Dwf, def.nu.Dws,...
    flip, TR, TEp, TEm,...
    'exchg', bool.exchg,...
    'mag', bool.mag);
  
  % 2-compartment dess signal gradients
  [spgradx, smgradx] = dess_2comp_gradx(...
    tmp.M0, tmp.ff, tmp.T1f, tmp.T1s, tmp.T2f, tmp.T2s,...
    def.nu.kap, def.nu.Dwf, def.nu.Dws,...
    flip, TR, TEp, TEm,...
    'exchg', bool.exchg,...
    'mag', bool.mag);
  
  % 2-compartment dess signal hessians
  [sphessx, smhessx] = dess_2comp_hessx(...
    tmp.M0, tmp.ff, tmp.T1f, tmp.T1s, tmp.T2f, tmp.T2s,...
    def.nu.kap, def.nu.Dwf, def.nu.Dws,...
    flip, TR, TEp, TEm,...
    'exchg', bool.exchg,...
    'mag', bool.mag);
  
  % plots
  figure(i);
  subplot(3,1,1);
  plot(tmp.(xf{i}), sp, tmp.(xf{i}), sm);
  ylabel('mag sig (a.u.)');
  legend('defocusing', 'refocusing');
  subplot(3,1,2);
  plot(tmp.(xf{i}), spgradx(:,i), tmp.(xf{i}), smgradx(:,i));
  ylabel(sprintf('1st deriv of mag sig w.r.t. %s (1/%s)', xf{i}, unit{i}));
  subplot(3,1,3);
  plot(tmp.(xf{i}), sphessx(:,i,i), tmp.(xf{i}), smhessx(:,i,i));
  ylabel(sprintf('2nd deriv of mag sig w.r.t. %s (1/%s^2)', xf{i}, unit{i}));
  xlabel(sprintf('%s (%s)', xf{i}, unit{i}));
end
