  function s = mese(m0, t1, t2, kap, ex, ref, ncyc, tr, te, varargin)
%|function s = mese(m0, t1, t2, kap, ex, ref, ncyc, tr, te, varargin)
%|
%|  mese signal model via extended phase graph simulation
%|    neglects incomplete transverse relaxation, so more accurate for longer tr
%|    neglects partial gradient spoiling due to bulk off-resonance effects
%|    neglects intrapulse relaxation, dephasing, etc.
%|    neglects interpulse exchange and diffusion
%|
%|  inputs
%|    m0          [(dim.o)]         spin density
%|    t1          [(dim.o)]         spin-lattice relaxation time            ms
%|    t2          [(dim.o)]         spin-spin relaxation time               ms
%|    kap         [(dim.o)]         flip angle scaling (same for ex/ref)
%|    ex
%|     .a         [1]               excitation nominal flip angle           rad
%|     .ph        [1]               excitation phase relative to +x axis    rad
%|    ref
%|     .a         [dim.e]           refocusing nominal flip angles          rad
%|     .ph        [dim.e]           refocusing phase relative to +x axis    rad
%|    ncyc        [dim.e 2]         across-voxel phase cycles 
%|                                    imparted by gradient crusher pairs
%|    tr          [1]               repetition time                         ms
%|    te          [dim.e]           echo intervals                          ms
%|  
%|  options
%|    mask        [(dim.o)]         object mask                             def: true(odims)
%|    mz0         [(dim.o)]         initial longitudinal magnetization      def: m0*(1-exp(-tr/t1))
%|    bool
%|     .mag       false|true        toggle magnitude signal off|on          def: false
%|     .video     false|true        toggle showing epg animation            def: false 
%|
%|  outputs
%|    s           [(dim.o) dim.e]   mese signals     
%|  
%|  copyright 2017, gopal nataraj, university of michigan
%|
%|  version control
%|    2017-10-30                    original
%|    2018-01-08                    added epg animation video            

% dimesions
dim.o = size(m0);
dim.e = length(te);

% default values
arg.mask = [];
arg.mz0 = [];
arg.bool.mag = 0;
arg.bool.video = 0;

% substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% if no mask specified, extrapolate to all voxels
if isempty(arg.mask)
  arg.mask = true(dim.o);
  dim.v = prod(dim.o);
else
  dim.v = numel(arg.mask(arg.mask));
end

% initial magnetization
if isempty(arg.mz0)
  tmp = bsxfun(@rdivide, tr, t1);
  tmp = exp(-tmp);
  tmp = bsxfun(@minus, 1, tmp);
  arg.mz0 = m0 .* tmp;
end

% vectorize 
m0 = masker(m0, arg.mask);
t1 = masker(t1, arg.mask);
t2 = masker(t2, arg.mask);
kap = masker(kap, arg.mask);
arg.mz0 = masker(arg.mz0, arg.mask);

% max number of nonzero configurations governed by sum of phase increments
tmp = sum(col(ncyc));
dim.c = 2*tmp+1;
config = zeros(3, dim.c, dim.v);

% config index + offset = array index
offset = tmp+1;

% initialize configurations
config(3, 0+offset, :) = arg.mz0;

% excitation pulse
config = rotx(ex.a, ex.ph, kap, config);

% iteratively evaluate echo amplitudes
s = nan(dim.v, dim.e);
for i = 1:dim.e
  % configurations relax separately for duration te(i)/2
  config = relax(m0, t1, t2, te(i)/2, config);
  
  % gradient lobe causes transverse configurations to increment
  config = precess(ncyc(i,1), config);
  
  % refocusing pulse
  config = rotx(ref.a(i), ref.ph(i), kap, config);
  
  % gradient lobe again causes transverse configurations to increment
  config = precess(ncyc(i,2), config);
  
  % configurations again relax separately for duration te(i)/2
  config = relax(m0, t1, t2, te(i)/2, config);
  
  % store ith echo amplitude
  s(:,i) = config(1, 0+offset, :);
  
  % video
  if arg.bool.video
    figure(1); 
    imagesc((1:dim.c)-offset, 1:3, abs(config));
    colorbar;
    caxis([0 0.8]);
    yticklabels([]);
    xlabel('configuration number', 'fontsize', 16);
    title(sprintf('echo %u: abs(a_0) = %0.4f', i, abs(s(:,i))));
    text(ones(3,1)*(-offset), col(1:3), {'abs(a_k)','abs(conj(a_k))','abs(b_k)'},...
      'VerticalAlignment', 'bottom',...
        'HorizontalAlignment', 'center',...
        'FontSize', 16,...
        'Color', 'k',...
        'Rotation', 90);
    pause(3/i);
  end
end

% option: take magnitude
if arg.bool.mag
  s = abs(s);
end

% embed s for output
s = embed(s, arg.mask);
end
