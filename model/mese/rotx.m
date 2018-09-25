  function out = rotx(ex, ph, kap, init)
%|function out = rotx(ex, ph, kap, init)
%|
%|  counterclockwise rotation operator, relative to +x-axis
%|
%|  inputs
%|    ex          [1]               nominal rotation angle                  rad
%|    ph          [1]               rotation phase relative to +x axis      rad
%|    kap         [dim.v]           rotation angle spatial variation
%|    init        [3 dim.c dim.v]   initial configurations
%|  
%|  output
%|    out         [3 dim.c dim.v]   rotated configurations
%|
%|  copyright 2017, gopal nataraj, university of michigan
%|
%|  version control
%|    2017-10-31                    original

% dimensions
dim.c = size(init,2);
dim.v = size(init,3);

% ccw rotation about +z axis by -ph*kap
tmp = exp(kap*(1i*ph));                             % [dim.v]
rotz = transpose([tmp 1./tmp ones(dim.v,1)]);       % [3 dim.v]
rotz = reshape(rotz, [3 1 dim.v]);                  % [3 1 dim.v]
out = bsxfun(@times, 1./rotz, init);                % [3 dim.c dim.v]

% ccw rotation about +x axis by +ex*kap
tmp = kap*ex;
c = cos(tmp);
s = sin(tmp);
tmp1 = (1+c)/2;
tmp2 = (1-c)/2;
tmp3 = s*(1i);
rotex = NaN(3,3,1,dim.v);
rotex(1,1,1,:) = tmp1;
rotex(1,2,1,:) = tmp2;
rotex(1,3,1,:) = -tmp3;
rotex(2,1,1,:) = tmp2;
rotex(2,2,1,:) = tmp1;
rotex(2,3,1,:) = tmp3;
rotex(3,1,1,:) = -tmp3/2;
rotex(3,2,1,:) = tmp3/2;
rotex(3,3,1,:) = c;
tmp4 = reshape(out, [1 3 dim.c dim.v]);
tmp4 = bsxfun(@times, rotex, tmp4);                 % [3 3 dim.c dim.v]
out = squeeze(sum(tmp4,2));                         % [3 dim.c dim.v]

% ccw rotation about +z axis +ph*kap
out = bsxfun(@times, rotz, out);                    % [3 dim.c dim.v]
end