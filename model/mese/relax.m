  function out = relax(m0, t1, t2, t, init)
%|function out = relax(m0, t1, t2, t, init)
%|
%|  relaxation operator
%|
%|  inputs
%|    m0          [dim.v]           spin density
%|    t1          [dim.v]           spin-lattice relaxation time            ms
%|    t2          [dim.v]           spin-spin relaxation time               ms
%|    t           [1]               relaxation duration                     ms
%|    init        [3 dim.c dim.v]   initial configuration         
%| 
%|  output
%|    out         [3 dim.c dim.v]   relaxed configuration
%|
%|  copyright 2017, gopal nataraj, university of michigan
%|
%|  version control
%|    2017-10-31                    original

% dimensions
dim.c = size(init,2);
dim.v = size(init,3);

% relaxation operation
e1 = exp(-t./t1);
e2 = exp(-t./t2);
tmp = transpose([e2 e2 e1]);                        % [3 dim.v]
tmp = reshape(tmp, [3, 1, dim.v]);                  % [3 1 dim.v]
out = bsxfun(@times, tmp, init);                    % [3 dim.c dim.v]

% recovery only affects coherent configuration
tmp = m0.*(1-e1);
tmp = reshape(tmp, [1 1 dim.v]);
out(3, (dim.c-1)/2+1, :) = out(3, (dim.c-1)/2+1, :) + tmp;
end