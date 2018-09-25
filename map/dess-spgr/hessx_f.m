  function [hess] = hessx_f(x, nu, P, dim, bool)
%|function [hess] = hessx_f(x, nu, P, dim, bool)
%|
%|  two-compartment signal model hessian evaluation w.r.t. latent object parameters x
%|
%|  inputs
%|    x         {L cell}        latent object parameters with cells size [V]
%|    nu        {N cell}        known parameters
%|    P         [1x1 struct]    scan parameters
%|                                1st field: data type (sp,de)
%|                                2nd field: scan parameter, as appropriate
%|     .*.tr    [S.*]               repetition times                                          ms
%|     .*.te    [S.* E.*]           echo times                                                ms
%|     .*.aex   [S.*]               nominal flip angle of excitation                          rad
%|    dim       [1x1 struct]    object containing dimension info  
%|    bool      [1x1 struct]    boolean variables
%|     .exchg   false|true        estimate exchange map                
%|     .mag.*   false|true        using magnitude (sp,de) data  
%|
%|  outputs
%|    hess      [DV LLV sparse] hessian evaluation
%|
%|  version control
%|    1.1       2016-10-26      original
%|    1.2       2016-10-27      project onto nonnegative orthant
%|    1.3       2018-01-22      removed hessian nonnegativity constraint
%|                              removed passing kfs as an optional argument

% initialize with tensor format for convenient indexing
dim.Vloc = numel(x{1});
hess = zeros(dim.D, dim.Vloc, dim.L, dim.L);                            % [D Vloc L L]
row = 1;

% 2-compartment spgr function calls
for s = 1:dim.S.sp
  hess(row,:,:,:) = spgr_2comp_hessx(...
    x{1:6}, nu{1:3},...
    P.sp.aex(s), P.sp.tr(s), P.sp.te(s,1),...
    'R2p_f', nu{4},...
    'R2p_s', nu{5},...
    'exchg', bool.exchg,...
    'mag', bool.mag.sp);
  row = row+1;
end

% 2-compartment dess function calls
for s = 1:dim.S.de
  [hess(row,:,:,:), hess(row+dim.S.de,:,:,:)] = dess_2comp_hessx(...
    x{1:6}, nu{1:3},...
    P.de.aex(s), P.de.tr(s), P.de.te(s,1), P.de.te(s,2),...
    'R2p_f', nu{4},...
    'R2p_s', nu{5},...
    'exchg', bool.exchg,...
    'mag', bool.mag.de);
  row = row+1;
end

% construct sparse block-diagonal hessian matrix
hess = reshape(hess, [dim.D dim.Vloc dim.L^2]);                         % [D Vloc L^2]
hess = permute(hess, [1 3 2]);                                          % [D L^2 Vloc]
hess = spmat(hess);                                                     % [DVloc L^2*Vloc sparse]
end