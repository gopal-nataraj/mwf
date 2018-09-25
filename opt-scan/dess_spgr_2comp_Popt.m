  function [P, fval, flag] = dess_spgr_2comp_Popt(P0, costArg, gradArg, varargin)
%|function [P, fval, flag] = dess_spgr_2comp_Popt(P0, costArg, gradArg, varargin)
%|
%|  parameter optimization for 2-compartment parameter estimation
%|    uses matlab's fmincon(...) for constrained zeroth-order optimization
%|    uses dess_spgr_2comp_cost(...) to evaluate cost function at iterates
%|    uses dess_spgr_2comp_costgrad(...) to evaluate cost function gradient
%|
%|    options available for optimizing w.r.t. (relative) standard deviation, (r)std
%|    assumes gaussian priors for latent object parameters
%|    assumes fs and ksf are fixed via constraints
%|    assumes kap (flip scaling), Dwf/s (off-res mean), R2pf/s (off-res bw) all known
%|    if exchg on/off, then L=7/6 latent parameters and K=5/6 known parameters
%|
%|  inputs
%|    P0        [1x1 struct]    initial parameter set
%|                                1st field: data type (sp,de)
%|                                2nd field: scan parameter, as appropriate
%|     .*.tr    [S.*]               repetition times                                          ms
%|     .*.aex   [S.*]               nominal flip angle of excitation                          rad
%|    costArg   {1 2*noptc}     option name-value pairs, passed to dess_spgr_2comp_cost()
%|    gradArg   {1 2*noptg}     option name-value pairs, passed to dess_spgr_2comp_costgrad()
%|
%|  options
%|    boxcon    [1x1 struct]    box constraints, [lower; upper]
%|                                1st arg: data type (sp,de)
%|                                2nd arg: scan parameter
%|     .sp.aex  [2]                 spgr flip angle                     def: [pi/180 pi/2]    rad
%|     .de.aex  [2]                 dess flip angle                     def: [pi/180 pi/2]    rad
%|     .sp.tr   [2]                 spgr repetition time                def: [11.8 Inf]       ms
%|     .de.tr   [2]                 dess repetition time                def: [17.5 Inf]       ms
%|    lincon    [1x1 struct]    linear inequality constraints
%|     .aex     [1]               total aex constraint                  def: Inf              rad
%|     .tr      [1]               total tr constraint                   def: sum(P0.*.tr)     ms
%|    fmincon   [1x1 struct]    fmincon options
%|     .wgrad   [char]            include cost function gradient        def: 'on'
%|     .disp    [char]            console feedback option               def: 'iter-detailed'
%|     .alg     [char]            optimization algorithm                def: 'active-set'
%|     .tolFun  [1]               cost function tolerance               def: 1e-8
%|     .tolX    [1]               iterate tolerance                     def: 1e-10
%|     .maxIter [1]               maximum number of iterations          def: 400
%|
%|  outputs
%|    P         [1x1 struct]    optimized parameter set (same structure as P0)
%|                                1st field: data type (sp,de)
%|                                2nd field: scan parameter, as appropriate
%|    fval      [1]             objective function value, at optimized P
%|    flag      [1]             exitflag from fmincon()
%|
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2016-04-08      original
%|    1.2       2016-04-16      added control over fmincon options
%|    2.1       2016-05-02      added option to include cost function gradient
%|    2.2       2016-08-17      changed format of P, grad, and boxcon

% constant declarations
S.de = length(P0.de.aex);
S.sp = length(P0.sp.aex);

% error checks on TRd, TRs
if length(P0.de.tr) ~= S.de
    error('flipd and TRd of unequal length!');
elseif length(P0.sp.tr) ~= S.sp
    error('flips and TRs of unequal length!');
end

% default box constraints
arg.boxcon.sp.aex   = col([1 90]) * pi/180;                             % rad
arg.boxcon.de.aex   = col([1 90]) * pi/180;                             % rad
arg.boxcon.sp.tr    = col([11.8 Inf]);                                  % ms
arg.boxcon.de.tr    = col([17.5 Inf]);                                  % ms

% default total linear constraints
arg.lincon.aex = Inf;                                                   % rad
arg.lincon.tr = Inf;                                                    % ms

% default optimization options
arg.fmincon.wgrad = 'on';     
arg.fmincon.disp = 'iter-detailed';
arg.fmincon.alg = 'interior-point';
arg.fmincon.tolFun = 1e-8;
arg.fmincon.tolX = 1e-10;
arg.fmincon.maxIter = 400;
  
% substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% convert P0 struct -> vec
P0v = [P0.de.aex; P0.sp.aex; P0.de.tr; P0.sp.tr];                       % [2(S.de+S.sp) 1]

% convert linear constraints
A = [ones(1,S.de+S.sp) zeros(1,S.de+S.sp);...
     zeros(1,S.de+S.sp) ones(1,S.de+S.sp)];                             % [2 2(S.de+S.sp)]
b = [arg.lincon.aex; arg.lincon.tr];                                    % [2 1]

% no equality constraints
Aeq = [];
beq = [];

% convert box constraints
lb = [...
  ones(S.de,1)*arg.boxcon.de.aex(1);...
  ones(S.sp,1)*arg.boxcon.sp.aex(1);...
  ones(S.de,1)*arg.boxcon.de.tr(1);...
  ones(S.sp,1)*arg.boxcon.sp.tr(1)];
ub = [...
  ones(S.de,1)*arg.boxcon.de.aex(2);...
  ones(S.sp,1)*arg.boxcon.sp.aex(2);...
  ones(S.de,1)*arg.boxcon.de.tr(2);...
  ones(S.sp,1)*arg.boxcon.sp.tr(2)];

% no nonlinear constraints
nonlcon = [];

% optimization options
opt = optimoptions('fmincon',...
    'GradObj', arg.fmincon.wgrad,...
    'Display', arg.fmincon.disp,...
    'Algorithm', arg.fmincon.alg,...
    'TolFun', arg.fmincon.tolFun,...
    'TolX', arg.fmincon.tolX,...
    'MaxIter', arg.fmincon.maxIter);

% anonymous function to handle extra arguments
fn = @(Pv) dess_spgr_2comp_Popt_helper(Pv, S, costArg, gradArg);

% constrained optimization
[Pv, fval, flag] = fmincon(fn, P0v, A, b, Aeq, beq, lb, ub, nonlcon, opt);  % [2(S.de+S.sp) 1]
        
% convert Pv vec -> P struct for output
P.de.aex  = Pv(1:S.de);
P.sp.aex  = Pv(S.de+1:S.de+S.sp);
P.de.tr   = Pv(S.de+S.sp+1:2*S.de+S.sp);
P.sp.tr   = Pv(2*S.de+S.sp+1:end);
end
  
  
  function [cost, gradv] = dess_spgr_2comp_Popt_helper(Pv, S, costArg, gradArg)
%|function [cost, gradv] = dess_spgr_2comp_Popt_helper(Pv, S, costArg, gradArg)
%| 
%|  helper function to dess_spgr_2comp_Popt
%|
%|  inputs
%|    Pv        [2(S.de+S.sp)]  vectorized parameters
%|    S         [1x1 struct]    scan number object
%|     .de      [1]               number of dess scans
%|     .sp      [1]               number of spgr scans
%|    costArg   {1 2*nopt}      option name-value pairs, passed to dess_spgr_2comp_cost()
%|    gradArg   {1 2*nopt}      option name-value pairs, passed to dess_spgr_2comp_costgrad()
%|
%|  outputs
%|    cost      [1]             output cost of dess_spgr_2comp_cost(...)
%|    gradv     [2(S.de+S.sp)]  vectorized output gradient of dess_spgr_2comp_costgrad(...)
%| 
%|  version control
%|    1.1       2016-04-08      original
%|    2.1       2016-05-02      added option to include cost function gradient
%|    2.2       2016-08-16      changed format of P and grad

% convert Pv vec -> P struct
P.de.aex = Pv(1:S.de,1);
P.sp.aex = Pv(S.de+1:S.de+S.sp);
P.de.tr   = Pv(S.de+S.sp+1:2*S.de+S.sp);
P.sp.tr   = Pv(2*S.de+S.sp+1:end);

% call cost function
cost = dess_spgr_2comp_cost(P, costArg{:});

% if required, call gradient function
if nargout > 1
  grad = dess_spgr_2comp_costgrad(P, gradArg{:});
  gradv = [grad.de.aex; grad.sp.aex; grad.de.tr; grad.sp.tr];               % [2(S.de+S.sp) 1]
end
end
