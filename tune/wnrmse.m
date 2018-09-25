  function out = wnrmse(xhat, x, varargin)
%|function out = wnrmse(xhat, x, varargin)
%|
%|  computes weighted normalized root mean squared latent parameter estimation error
%|
%|  inputs
%|    xhat      [1x1 struct]    latent parameter estimate
%|     .m0      [(odims)]         spin density
%|     .ff      [(odims)]         fast fraction
%|     .t1f     [(odims)]         fast spin-lattice relaxation time                           ms
%|     .t1s     [(odims)]         slow spin-lattice relaxation time                           ms
%|     .t2f     [(odims)]         fast spin-spin relaxation time                              ms
%|     .t2s     [(odims)]         slow spin-spin relaxation time                              ms
%|    x         [1x1 struct]    latent parameter reference
%|     .m0      [(odims)]         spin density
%|     .ff      [(odims)]         fast fraction
%|     .t1f     [(odims)]         fast spin-lattice relaxation time                           ms
%|     .t1s     [(odims)]         slow spin-lattice relaxation time                           ms
%|     .t2f     [(odims)]         fast spin-spin relaxation time                              ms
%|     .t2s     [(odims)]         slow spin-spin relaxation time                              ms
%|
%|  options
%|    mask      [(odims)]       defines roi on which to compute error   def: true(odims)
%|    wght      [L]             relative weighting of parameters        def: ones(L)/L
%|
%|  output
%|    out       [1]             weighted normalized root mean squared error
%|
%|  copyright 2017, gopal nataraj, university of michigan
%|  
%|  version control
%|    1.1       2018-02-23      adapted from m0,t1,t2 estimation version

% dimensions
dim.odims = size(xhat.m0);
dim.L = length(fieldnames(xhat));

% default values
arg.mask = true(dim.odims);
arg.wght = ones(dim.L,1);

% substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% convert structs to cell arrays
xhat = struct2cell(xhat);
x = struct2cell(x);

% vectorize inputs
for l = 1:dim.L
  xhat{l} = masker(xhat{l}, arg.mask);                                    % [V]
  x{l} = masker(x{l}, arg.mask);                                          % [V]
end
dim.V = length(xhat{1});

% normalize weights
tmp = div0(1,sum(arg.wght));
if abs(tmp-1)>eps
  warn('weights should add to 1; normalizing...');
  arg.wght = arg.wght * tmp;
end
  
% weighted nrmse
num = abs([xhat{:}]-[x{1:dim.L}]);                                        % [V L]
den = bsxfun(@max, abs([x{1:dim.L}]), eps);                               % [V L]
tmp = transpose(div0(num,den));                                           % [L V]
tmp = tmp.^2;                                                             % [L V]
tmp = diag(arg.wght) * tmp;                                               % [L V]
out = col(sum(tmp, 1));                                                   % [V]
out = mean(out);                                                          % [1]
out = sqrt(out);                                                          % [1]
end
