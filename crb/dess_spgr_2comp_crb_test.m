% script dess_spgr_2comp_crb_test.m
% header file to compute fisher info from dess/spgr combinations
%
% copyright 2016, gopal nataraj, university of michigan
%
% version control
%   1.1     2016-03-31      original
%   1.2     2016-04-08      add backslash, trace to comp w/ dess_spgr_2comp_cost.m
%   1.3     2016-04-15      modified to test fast no-exchange mode

% setup
if (~exist('irtdir', 'var'))
  curdir = cd('~/Box Sync/work/irt');
  irtdir = pwd;
  setup(); 
  cd(curdir);
end
addpath('../model/spgr/');
addpath('../model/dess/');

% unconstrained, unknown object parameters
M0 = 1;
ff = 0.20;
T1f = 400;                                          % ms
T1s = 1000;                                         % ms
T2f = 20;                                           % ms
T2s = 80;                                           % ms

% known object parameters
kap = 1;
Dwf = 0;                                            % kHz
Dws = 0;                                            % kHz

% noise standard deviation
oldVoxelVol = 1 * 1 * 5;                            % mm^3
newVoxelVol = (240/256) * (240/256) * (30/6);       % mm^3
oldvar_im_1coil = 2.62e-7;                          % a.u.; at 1x1x5 resolution 
newvar_im_1coil = oldvar_im_1coil * (oldVoxelVol / newVoxelVol);
noise_var_ssos = newvar_im_1coil/2;                 % hi-snr approx
noise_std_ssos = sqrt(noise_var_ssos);

% dess scan parameters
nfd = 10;
nTRd = 1;
flipd_deg = linspace(1,180,nfd)';                   % deg
flipd = col(repmat((pi/180)*flipd_deg, [1 nTRd]));  % rad
TRd = col(repmat(linspace(15,40,nTRd), [nfd 1]));   % ms
TEp = 4.67 * col(ones(nfd,nTRd));                   % ms
TEm = 4.67 * col(ones(nfd,nTRd));                   % ms
sigp = noise_std_ssos * col(ones(nfd,nTRd));
sigm = noise_std_ssos * col(ones(nfd,nTRd));

% spgr scan parameters
nfs = 10;
nTRs = 1;
flips_deg = linspace(1,180,nfs)';                   % deg
flips = col(repmat((pi/180)*flips_deg, [1 nTRs]));  % rad
TRs = col(repmat(linspace(10,30,nTRs), [nfs 1]));   % ms
TEs = 4.67 * col(ones(nfs,nTRs));                   % ms
sigs = noise_std_ssos * col(ones(nfs,nTRs));        

% weighting matrices
W7 = diag([0; 1; 0; 0; 0; 0; 0]);
W6 = diag([0; 1; 0; 0; 0; 0]);

% output parameters
pr = 0;

% fisher information from all scans, neglecting exchange
tic;
F6d = dess_2comp_crb(ff, T1f, T1s, T2f, T2s, kap, Dwf, Dws,...
    flipd, TRd, TEp, TEm, sigp, sigm, 'exchg', false);
F6s = spgr_2comp_crb(ff, T1f, T1s, T2f, T2s, kap, Dwf, Dws,...
    flips, TRs, TEs, sigs, 'exchg', false);
t_exchg0 = toc;
F6 = squeeze(F6d + F6s);
varff6 = trace(W6 * (F6 \ (W6')));
sigff6 = sqrt(varff6);
cov6 = inv(F6);

% display
fprintf('Found unity-M0 ff relative std dev = %0.5f, in %0.3f ms.\n',...
  div0(sigff6,ff), t_exchg0*1000);
figure; im('notick', log10(abs(inv(squeeze(F6d)))), [-5 3], 'cbar', ' ');
figure; im('notick', log10(abs(inv(squeeze(F6s)))), [-5 3], 'cbar', ' ');
figure; im('notick', log10(abs(cov6)), [-5 3], 'cbar', ' ');
if (pr), print('-deps', 'dess_spgr_2comp_cov6x6.eps'), end;
