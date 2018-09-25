% script dess_2comp_sym.m
% symbolic toolbox calculation of two-compartment dess analytical model
%
% copyright 2016, university of michigan
%   1.1     2016-03-18      original
%   1.2     2016-03-29      split xy/z matrix exponentials
%   1.3     2016-03-30      added gradient, flip angle spatial variation
%   1.4     2016-04-14      added option to remove exchange, mixed partials
%   1.5     2016-06-15      changed sign of second echo; then changed back
%   1.6     2016-10-26      added hessian evaluation
%   2.1     2018-01-15      now separately considers real vs. complex cases
%   2.2     2018-01-16      replaced assume() calls with assumeAlso() calls
%   2.3     2018-01-17      now uses simpler form for mag sig with Dwf==Dws
%   2.4     2018-01-19      symbolic sign error workaround: manually take signal mag
%   2.5     2018-01-21      now defines dependent variable fs := 1-ff

% options
bool.exchg = 0;             % model exchange
bool.mag = 0;               % model magnitude signal and its derivatives
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
syms flip TR TEp TEm t positive;
syms phif phis real;
syms R2pf R2ps positive;

% assumptions
assumeAlso(ff <= 1);
assumeAlso(0 < kap & kap < 1);
assumeAlso(0 < flip & flip < pi/2);
assumeAlso(0 < TR);
assumeAlso(0 < TEp & TEp < TR);
assumeAlso(0 < TEm & TEm < TR);
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

% assemble and exponentiate transverse bloch matrix
if bool.exchg
  fprintf('Computing transverse matrix exponential, with exchange...');
  Axy =  [-1/T2f-kfs,     ksf,            Dwf,            0;      
          kfs,            -1/T2s-ksf,     0,              Dws;
          -Dwf,           0,              -1/T2f-kfs,     ksf;
          0,              -Dws,           kfs,            -1/T2s-ksf];
else
  fprintf('Computing transverse matrix exponential, without exchange...');
  Axy =  [-1/T2f,         0,              Dwf,            0;      
          0,              -1/T2s,         0,              Dws;
          -Dwf,           0,              -1/T2f,         0;
          0,              -Dws,           0,              -1/T2s];
end
[Vxy(t), Dxy(t)] = eig(Axy * t);
expAxy(t) = Vxy(t) * diag(exp(diag(Dxy(t)))) / Vxy(t);
fprintf('done.\n');

% assemble and exponentiate longitudinal bloch matrix
if bool.exchg
  fprintf('Computing longitudinal matrix exponential, with exchange...');
  Az =   [-1/T1f-kfs,     ksf;
          kfs,            -1/T1s-ksf];
else
  fprintf('Computing longitudinal matrix exponential, without exchange...');
  Az =   [-1/T1f,         0;
          0,              -1/T1s];
end
[Vz(t), Dz(t)] = eig(Az * t);
expAz(t) = Vz(t) * diag(exp(diag(Dz(t)))) / Vz(t);
fprintf('done.\n');

% construct matrix exponential blockwise
expA(t) =  [expAxy(t),  zeros(4,2);...
            zeros(2,4), expAz(t)];

% assemble constant matrix
c = M0 * [zeros(4,1); ff/T1f; fs/T1s];

% assemble rotation matrix
R6 = [1,      0,      0,      0,      0,      0;
      0,      1,      0,      0,      0,      0;
      0,      0,      cos(a), 0,      sin(a)  0;
      0,      0,      0,      cos(a), 0       sin(a);
      0,      0,      -sin(a),0,      cos(a), 0;
      0,      0,      0,      -sin(a),0,      cos(a)];

% steady-state magnetization, just after rf pulse
A =    [Axy, zeros(4,2);
        zeros(2,4), Az];
mp_0 = ((eye(6) - R6*expA(TR)) \ R6) * (expA(TR) - eye(6)) * (A \ c);
mm_0 = R6 \ mp_0;

% with exchange, integrate transverse signals
if bool.exchg
  return;
  
  % compartmental transverse magnetizations just after rf pulse
  simp_arg = {...
    'IgnoreAnalyticConstraints', true,...
    'Criterion', 'preferReal',...
    'Steps', 100};
  mpxy_0_f = simplify(...
    subs(mp_0(1) + 1i*mp_0(3), {Dwf, Dws}, {phif/TR, phis/TR}), simp_arg{:});
  mpxy_0_s = simplify(...
    subs(mp_0(2) + 1i*mp_0(4), {Dwf, Dws}, {phif/TR, phis/TR}), simp_arg{:});
  mmxy_0_f = simplify(...
    subs(mm_0(1) + 1i*mm_0(3), {Dwf, Dws}, {phif/TR, phis/TR}), simp_arg{:});
  mmxy_0_s = simplify(...
    subs(mm_0(2) + 1i*mm_0(4), {Dwf, Dws}, {phif/TR, phis/TR}), simp_arg{:});

  % compartmental transverse signals just after rf pulse
  scale = div0(1, (2*pi)^2);
  int_arg = {...
    'IgnoreAnalyticConstraints', true,...
    'IgnoreSpecialCases', true,...
    'PrincipalValue', true...
  };
  spxy_0_f = scale * int(int(mpxy_0_f, phif, 0, 2*pi, int_arg{:}), phis, 0, 2*pi, int_arg{:});
  spxy_0_s = scale * int(int(mpxy_0_s, phif, 0, 2*pi, int_arg{:}), phis, 0, 2*pi, int_arg{:});
  smxy_0_f = scale * int(int(mmxy_0_f, phif, 0, 2*pi, int_arg{:}), phis, 0, 2*pi, int_arg{:});
  smxy_0_s = scale * int(int(mmxy_0_s, phif, 0, 2*pi, int_arg{:}), phis, 0, 2*pi, int_arg{:});
% without exchange, use explicit expressions
else
  E1f = exp(-TR / T1f);
  v1f = (1 - E1f*cos(a)) / (E1f - cos(a));
  E2f = exp(-TR / T2f);
  Rf = sqrt((1-E2f^2) / (1-E2f^2/v1f^2));
  spxy_0_f = +1i*M0*ff*tan(a/2) * (1 - Rf/v1f);
  smxy_0_f = -1i*M0*ff*tan(a/2) * (1 - Rf);

  E1s = exp(-TR / T1s);
  v1s = (1 - E1s*cos(a)) / (E1s - cos(a));
  E2s = exp(-TR / T2s);
  Rs = sqrt((1-E2s^2) / (1-E2s^2/v1s^2));
  spxy_0_s = +1i*M0*fs*tan(a/2) * (1 - Rs/v1s);
  smxy_0_s = -1i*M0*fs*tan(a/2) * (1 - Rs);
end

% compartment-wise dephasing between rf and signal reception
spxy_tep_f = spxy_0_f * exp(-TEp * (R2pf + 1/T2f)) * exp(+1i * TEp * Dwf);
spxy_tep_s = spxy_0_s * exp(-TEp * (R2ps + 1/T2s)) * exp(+1i * TEp * Dws);
smxy_tem_f = smxy_0_f * exp(-TEm * (R2pf - 1/T2f)) * exp(-1i * TEm * Dwf);
smxy_tem_s = smxy_0_s * exp(-TEm * (R2ps - 1/T2s)) * exp(-1i * TEm * Dws);

% total received signals
if bool.mag && ~bool.freq
  spxy_tep = abs(spxy_tep_f) + abs(spxy_tep_s);
  smxy_tem = abs(smxy_tem_f) + abs(smxy_tem_s);
elseif bool.mag
  spxy_tep = abs(spxy_tep_f + spxy_tep_s);
  smxy_tem = abs(smxy_tem_f + smxy_tem_s);
else
  spxy_tep = spxy_tep_f + spxy_tep_s;
  smxy_tem = smxy_tem_f + smxy_tem_s;
end

% function suffix
suffix = sprintf('_exchg%u_mag%u_freq%u', bool.exchg, bool.mag, bool.freq);

% variables
if bool.mag && ~bool.freq
  var1 = [M0 ff T1f T1s T2f T2s kap flip TR TEp R2pf R2ps];
  var2 = [M0 ff T1f T1s T2f T2s kap flip TR TEm R2pf R2ps];
else
  var1 = [M0 ff T1f T1s T2f T2s kap Dwf Dws flip TR TEp R2pf R2ps];
  var2 = [M0 ff T1f T1s T2f T2s kap Dwf Dws flip TR TEm R2pf R2ps];
end

% signal models
tmp = ['dess_2comp_echo1', suffix, '.m'];
if ~exist(tmp, 'file')
  fpxy_tep = matlabFunction(spxy_tep,...
    'file', tmp,...
    'vars', var1);
end
tmp = ['dess_2comp_echo2', suffix, '.m'];
if ~exist(tmp, 'file')
  fmxy_tem = matlabFunction(smxy_tem,...
    'file', tmp,...
    'vars', var2);
end

% row gradients of signal models w.r.t. x
tmp = ['dess_2comp_echo1', suffix, '_gradx.m'];
if ~exist(tmp, 'file')
  fprintf('Computing row gradient of first dess signal w.r.t. x...');
  spxy_tep_gradx = simplify(jacobian(spxy_tep, [M0 ff T1f T1s T2f T2s]));   % [1 L]
  fpxy_tep_gradx = matlabFunction(spxy_tep_gradx,...
    'file', tmp,...
    'vars', var1);
  fprintf('done.\n');
end
tmp = ['dess_2comp_echo2', suffix, '_gradx.m'];
if ~exist(tmp, 'file')
  fprintf('Computing row gradient of second dess signal w.r.t. x...');
  smxy_tem_gradx = simplify(jacobian(smxy_tem, [M0 ff T1f T1s T2f T2s]));   % [1 L]
  fmxy_tem_gradx = matlabFunction(smxy_tem_gradx,...
    'file', tmp,...
    'vars', var2);
  fprintf('done.\n');
end

% mixed gradients of signal models w.r.t. x, P
tmp = ['dess_2comp_echo1', suffix, '_gradx_gradp.m'];
if ~exist(tmp, 'file')
  if ~exist('spxy_tep_gradx', 'var')
    spxy_tep_gradx = simplify(jacobian(spxy_tep, [M0 ff T1f T1s T2f T2s])); % [1 L]
  end
  fprintf('Computing mixed gradient of first dess signal w.r.t. x, P...');
  spxy_tep_gradx = transpose(spxy_tep_gradx);                               % [L]
  spxy_tep_gradx_gradp = simplify(jacobian(spxy_tep_gradx, [flip TR]));     % [L P]
  fpxy_tep_gradx_gradp = matlabFunction(spxy_tep_gradx_gradp,...
    'file', tmp,...
    'vars', var1);
  fprintf('done.\n');
end
tmp = ['dess_2comp_echo2', suffix, '_gradx_gradp.m'];
if ~exist(tmp, 'file')
  if ~exist('smxy_tem_gradx', 'var')
    smxy_tem_gradx = simplify(jacobian(smxy_tem, [M0 ff T1f T1s T2f T2s])); % [1 L]
  end
  fprintf('Computing mixed gradient of second dess signal w.r.t. x, P...');
  smxy_tem_gradx = transpose(smxy_tem_gradx);                               % [L]
  smxy_tem_gradx_gradp = simplify(jacobian(smxy_tem_gradx, [flip TR]));     % [L P]
  fmxy_tem_gradx_gradp = matlabFunction(smxy_tem_gradx_gradp,...
    'file', tmp,...
    'vars', var2);
  fprintf('done.\n');
end
  
% hessians of signal models w.r.t. x
tmp = ['dess_2comp_echo1', suffix, '_hessx.m'];
if ~exist(tmp, 'file')
  fprintf('Computing hessian of first dess signal w.r.t. x...');
  spxy_tep_hessx = simplify(hessian(spxy_tep, [M0 ff T1f T1s T2f T2s]));    % [L L]
  spxy_tep_hessx = reshape(spxy_tep_hessx, [1 6 6]);                        % [1 L L]
  fpxy_tep_hessx = matlabFunction(spxy_tep_hessx,...
    'file', tmp,...
    'vars', var1);
  fprintf('done.\n');
end
tmp = ['dess_2comp_echo2', suffix, '_hessx.m'];
if ~exist(tmp, 'file')
  fprintf('Computing hessian of second dess signal w.r.t. x...');
  smxy_tem_hessx = simplify(hessian(smxy_tem, [M0 ff T1f T1s T2f T2s]));    % [L L]
  smxy_tem_hessx = reshape(smxy_tem_hessx, [1 6 6]);                        % [1 L L]
  fmxy_tem_hessx = matlabFunction(smxy_tem_hessx,...
    'file', tmp,...
    'vars', var2);
  fprintf('done.\n');
end
