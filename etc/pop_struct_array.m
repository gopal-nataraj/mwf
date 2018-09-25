  function array = pop_struct_array(fn, elem, dims, arg)
%|function array = pop_struct_array(fn, elem, dims, arg)
%|
%|  evaluates a function element-wise on a struct array
%|    requires recursive calls to pop_struct_array_helper(...)
%|
%|  input
%|    fn        [func handle]   element-wise function operator
%|                                assumes form fn(idx, arg) 
%|                                evaluates on array at row vector position idx
%|    elem      [struct]        sample array element w/ (possibly empty) fields defined
%|    dims      [1 ndims]       output array dimensions
%|    arg       [struct]        obj containing add'l params needed in fn
%|
%|  output
%|    array     [(dims)]        struct array, w/ fn evaluated at every possible idx
%|
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2016-04-20      modified populate_array(...)

% instantiate struct array
array = repmat(elem, dims);     % [(odims) struct]

% root call to recursive helper
array = pop_struct_array_helper(fn, array, dims, dims, [], 1, arg);
end


  function array = pop_struct_array_helper(fn, array, subdims, fulldims, idx, curr, arg)
%|function array = pop_struct_array_helper(fn, array, subdims, fulldims, idx, curr, arg)
%|
%|  helper function to pop_struct_array(...)
%|
%|  input
%|    fn        [func handle]   element-wise function operator
%|                                assumes form fn(idx, arg) 
%|                                evaluates on array at row vector position idx
%|    array     [(dims)]        output array, progressively populated
%|    subdims   [1 0-ndims]     subset of output array dimensions
%|    fulldims  [1 ndims]       full output array dimensions
%|    idx       [1 0-ndims]     subset of current index at which to eval fn
%|    curr      [1]             which dim currently being iterated
%|    arg       [struct]        obj containing add'l params needed in fn
%|
%|  version control
%|    1.1       2016-04-20      modified populate_array_helper(...)

% base case: evaluate fn, with idx of size [1 ndims]
if isempty(subdims)
  idx_c = num2cell(idx);        % {1 L+K}
  array(sub2ind(fulldims, idx_c{:})) = fn(idx, arg);
% recursive step: continue to build idx of size [1 <ndims]
else
  rest = subdims(2:end);
  for d = 1:subdims(1)
    if length(idx) < curr
      idx = [idx d];
    else
      idx(curr) = d;
    end
    array = pop_struct_array_helper(fn, array, rest, fulldims, idx, curr+1, arg);
  end
end
end