% script spgr_2comp_sym.m
% two-compartment spgr analytical model with exchange
%
% copyright 2016, university of michigan
%   1.1     2016-03-29      original
%   1.2     2016-03-30      added gradient, flip angle spatial variation
%   1.3     2016-04-15      added option to remove exchange, mixed partials
%   1.4     2016-10-26      added hessian evaluation
%   2.1     2018-01-21      incorporated v2.* changes from ../dess/dess_2comp_sym.m

% options
bool.exchg = 0;             % model exchange
bool.mag = 1;               % model magnitude signal and its derivatives
bool.freq = 0;              % model different compartmental off-resonance frequencies 

% symbolic variable declarations
syms M0 ff positive;
syms T1f T1s T2f T2s positive;
if bool.exchg
  syms kfs positive;
end
syms kap positive;
if bool.mag && ~bool.freq
  Dwf = 0;
  Dws = 0;
else
  syms Dwf Dws real;
end
syms flip TR TE t positive;
syms R2pf R2ps positive;

% assumptions
assumeAlso(ff <= 1);
assumeAlso(0 < kap & kap < 1);
assumeAlso(0 < flip & flip < pi/2);
assumeAlso(0 < TR);
assumeAlso(0 < TE & TE < TR);
assumeAlso(0 < t & t < TR);
assumeAlso(exp(-TR/T1f) ~= cos(flip*kap));
assumeAlso(exp(-TR/T1s) ~= cos(flip*kap));
if ~bool.mag || bool.freq
  assumeAlso(Dwf == Dws);
end

% slow-relaxing fraction
fs = 1-ff;

% slow-to-fast exchange rate
if bool.exchg
  ksf = kfs*(ff/fs);
end

% add spatial variation to flip
a = flip*kap;

% assemble and exponentiate longitudinal bloch matrix
if bool.exchg
  fprintf('Computing longitudinal matrix exponential, with exchange...');
  Az =  [-1/T1f-kfs,     ksf;
         kfs,            -1/T1s-ksf];
else
  fprintf('Computing longitudinal matrix exponential, without exchange...');
  Az =  [-1/T1f,         0;
         0,              -1/T1s];
end
[Vz(t), Dz(t)] = eig(Az * t);
expAz(t) = Vz(t) * diag(exp(diag(Dz(t)))) / Vz(t);
fprintf('done!\n');

% spoiled, full matrix
S = [zeros(4,4), zeros(4,2);...
     zeros(2,4), eye(2)];
SexpA(t) = [zeros(4,4), zeros(4,2);...
            zeros(2,4), expAz(t)];
        
% assemble constant vector
invAc = [zeros(4,1);...
         Az \ (M0 * [ff/T1f; fs/T1s])];

% assemble rotation matrix
R6 = [1,      0,      0,      0,      0,      0;
      0,      1,      0,      0,      0,      0;
      0,      0,      cos(a), 0,      sin(a)  0;
      0,      0,      0,      cos(a), 0       sin(a);
      0,      0,      -sin(a),0,      cos(a), 0;
      0,      0,      0,      -sin(a),0,      cos(a)];

% steady-state magnetization, just after rf pulse
ms_0 = ((eye(6) - R6*SexpA(TR)) \ R6) * (SexpA(TR) - S) * invAc;

% compartmental transverse signals just after rf pulse
ssxy_0_f = ms_0(1) + 1i*ms_0(3);
ssxy_0_s = ms_0(2) + 1i*ms_0(4);

% apply compartment-wise dephasing between tip and signal collection
ssxy_te_f = ssxy_0_f * exp(-TE * (R2pf + 1/T2f)) * exp(1i * TE * Dwf);
ssxy_te_s = ssxy_0_s * exp(-TE * (R2ps + 1/T2s)) * exp(1i * TE * Dws);

% total received signal
if bool.mag && ~bool.freq
  ssxy_te = abs(ssxy_te_f) + abs(ssxy_te_s);
elseif bool.mag
  ssxy_te = abs(ssxy_te_f + ssxy_te_s);
else
  ssxy_te = ssxy_te_f + ssxy_te_s;
end

% function suffix
suffix = sprintf('_exchg%u_mag%u_freq%u', bool.exchg, bool.mag, bool.freq);

% latent parameters
if bool.exchg
  x = [M0 ff T1f T1s T2f T2s kap];
else
  x = [M0 ff T1f T1s T2f T2s];
end
L = length(x);

% input variables
if bool.mag && ~bool.freq
  var = [M0 ff T1f T1s T2f T2s kap flip TR TE R2pf R2ps];
else
  var = [M0 ff T1f T1s T2f T2s kap Dwf Dws flip TR TE R2pf R2ps];
end

% signal model
tmp = ['spgr_2comp', suffix, '.m'];
if ~exist(tmp, 'file')
  fsxy_te = matlabFunction(ssxy_te,...
    'file', tmp,...
    'vars', var);
end

% row gradient of signal model w.r.t. x
tmp = ['spgr_2comp', suffix, '_gradx.m'];
if ~exist(tmp, 'file')
  fprintf('Computing row gradient of spgr signal w.r.t. x...');
  ssxy_te_gradx = simplify(jacobian(ssxy_te, x));                               % [1 L]
  fsxy_te_gradx = matlabFunction(ssxy_te_gradx,...
    'file', tmp,...
    'vars', var);
  fprintf('done.\n');
end

% mixed gradient of signal model w.r.t. x, P
tmp = ['spgr_2comp', suffix, '_gradx_gradp.m'];
if ~exist(tmp, 'file')
  if ~exist('ssxy_te_gradx', 'var')
    ssxy_te_gradx = simplify(jacobian(ssxy_te, x));                             % [1 L]
  end
  fprintf('Computing mixed gradient of spgr signal w.r.t x, P...');
  ssxy_te_gradx = transpose(ssxy_te_gradx);                                     % [L]
  ssxy_te_gradx_gradp = simplify(jacobian(ssxy_te_gradx, [flip TR]));           % [L P]
  fsxy_te_gradx_gradp = matlabFunction(ssxy_te_gradx_gradp,...
    'file', tmp,...
    'vars', var);
  fprintf('done.\n');
end

% hessian of signal model w.r.t. x, P
tmp = ['spgr_2comp', suffix, '_hessx.m'];
if ~exist(tmp, 'file')
  fprintf('Computing hessian of spgr signal w.r.t. x...');
  ssxy_te_hessx = simplify(hessian(ssxy_te, x));                                % [L L]
  ssxy_te_hessx = reshape(ssxy_te_hessx, [1 L L]);                              % [1 L L]
  fsxy_te_hessx = matlabFunction(ssxy_te_hessx,...
    'file', tmp,...
    'vars', var);
  fprintf('done.\n');
end
