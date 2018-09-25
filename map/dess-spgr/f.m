  function [fval] = f(x, nu, P, dim, bool)
%|function [fval] = f(x, nu, P, dim, bool)
%|
%|  two-compartment signal model function evaluation
%|
%|  inputs
%|    x         {L cell}        latent object parameters with cells size [V]
%|    nu        {N cell}        known parameters with cells size [Vloc]
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
%|    fval      [DV]            signal model function evaluations
%|
%|  version control
%|    1.1       2016-08-15      adapted from mri_m0t1t2inveff_map(...)
%|    1.2       2018-01-22      removed passing kfs as an optional argument

% initialize with matrix format for convenient indexing
dim.Vloc = numel(x{1});                                                
fval = NaN(dim.D, dim.Vloc);
row = 1;

% 2-compartment spgr function calls
for s = 1:dim.S.sp
  fval(row,:) = spgr_2comp(...
    x{1:6}, nu{1:3},...
    P.sp.aex(s), P.sp.tr(s), P.sp.te(s,1),...
    'R2p_f', nu{4},...
    'R2p_s', nu{5},...
    'exchg', bool.exchg,...
    'mag', bool.mag.sp);
  row = row+1;
end

% 2-compartment dess function calls
for s= 1:dim.S.de
  [fval(row,:), fval(row+dim.S.de,:)] = dess_2comp(...
    x{1:6}, nu{1:3},...
    P.de.aex(s), P.de.tr(s), P.de.te(s,1), P.de.te(s,2),...
    'R2p_f', nu{4},...
    'R2p_s', nu{5},...
    'exchg', bool.exchg,...
    'mag', bool.mag.de);
  row = row+1;
end

% reshape for output
fval = col(fval);                                                       % [DVloc]
end