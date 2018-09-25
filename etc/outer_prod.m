  function prod = outer_prod(vecs)
%|function prod = outer_prod(vecs)
%|
%|  evaluates an outer product over a series of vectors
%|
%|  input
%|    vecs      {N}             column vectors of length D_1...D_N
%|
%|  output
%|    prod      [(D_1...D_N)]   outer product
%|
%|  version control
%|    1.1       2016-04-07      original

% constant declarations
N = length(vecs);
prod = vecs{1};
D = length(vecs{1});

% construct outer product as a column first
for n = 2:N
  prod = kron(vecs{n}, prod);   % [prod(D_1,...,D_n)]
  D = [D length(vecs{n})];      % [1 n]
end

% reshape for output
prod = reshape(prod, D);
end