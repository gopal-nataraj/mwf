  function [Xkk_inv] = diag_inv(X, k, varargin)
%|function [Xkk_inv] = diag_inv(X, k, varargin)
%|
%|  efficiently computes diagonal elements of the inverse of a symmetric matrix
%|    uses schur complement formula
%|
%|  inputs
%|    X         [n n]       symmetric matrix
%|    k         [1]         diagonal element of inverse desired
%|
%|  options
%|    symcheck  [1]         threshold by which to check symmetry                def: 100*eps
%|
%|  outputs:
%|    Xkk_inv   [1]         inv(X)_{k,k}
%|
%|  copyright 2014, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2014-07-15  original
%|    1.2       2016-04-03  error checks, div0, backslash

% default values
arg.symcheck = 100*eps;

% substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% make sure X is symmetric
if norm(X - X',2) > arg.symcheck
    error('X is not symmetric!');
end

% make sure k is within range
if k<1 || k>size(X,1)
    error('k out of range!');
end

% extract submatrix and subvector
X_k =  [X(1:k-1,1:k-1)      X(1:k-1,k+1:end);...
        X(k+1:end,1:k-1)    X(k+1:end,k+1:end)];
a_k =  [X(1:k-1,k);... 
        X(k+1:end,k)];

% submatrix inversion
tmp = a_k' * (X_k \ a_k);
Xkk_inv = div0(1, X(k,k) - tmp);
  end