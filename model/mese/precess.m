  function out = precess(ncyc, init)
%|function out = precess(ncyc, init)
%|
%|  propagates epg configurations through an integral number of phase cycles
%|
%|  inputs
%|    ncyc        [1]               number of phase cycles 
%|    init        [3 dim.c dim.v]   initial configurations
%|
%|  outputs
%|    out         [3 dim.c dim.v]   updated configuration
%|
%|  copyright 2017, gopal nataraj, university of michigan
%|
%|  version control
%|    2017-10-31                    original

% dimensions
dim.c = size(init,2);
dim.v = size(init,3);

% check to ensure full phase cycles
if mod(ncyc,1)~=0
  error('ncyc must be an integer.');
end

% increment transverse states only
out = nan(3, dim.c, dim.v);
out(1,:,:) = circshift(init(1,:,:), [0 +ncyc 0]);
out(2,:,:) = circshift(init(2,:,:), [0 -ncyc 0]);
out(3,:,:) = init(3,:,:);
end