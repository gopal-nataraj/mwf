  function newCell = cell_rep(oldCell, varargin)
%|function newCell = cell_rep(oldCell, varargin)
%|  
%|  replaces {'name', val} pairs with new values
%|    useful for passing through selectively-modified variable arguments to subfunctions,
%|
%|  input
%|    oldCell   {1 2*nopts}     original {'name', value} pairs
%|
%|  output
%|    newCell   {1 >=2*nopts}   new {'name', value} pairs, possibly with additions
%| 
%|  copyright 2016, gopal nataraj, university of michigan
%|  
%|  version control
%|    1.1       2016-04-20      original

% check to make sure varargin has paired inputs
if rem(length(varargin),2) ~= 0
  error('Need {''name'', val} pairs as input.');
end

% instantiation
newCell = oldCell;
names = varargin(1:2:end);
vals = varargin(2:2:end);

% replace first instance of 'name' with corresponding new value
for i = 1:length(names)
  for j = 1:2:length(oldCell)
    if streq(newCell{j}, names{i})
      newCell{j+1} = vals{i};
      break;
    elseif j == length(oldCell)-1
      warning('variable %s not found. Adding to end of newCell.', names{i});
      newCell = {newCell{:}, names{i}, vals{i}};
    end
  end
end

end
  