% script mese_test.m
% test script for mese.m
%
% gopal nataraj
% copyright 2017, university of michigan
%
% version control
%   2017-10-31      original

% irt setup
if ~exist('irtdir', 'var')
  curdir = cd('~/Box Sync/work/irt');
  irtdir = pwd;
  setup();
  cd(curdir);
end

% perfect excitation, imperfect refocusing
m0 = 1;
t1 = 800;
t2 = 50;
kap = 1;

ex.a = pi/2;
ex.ph = pi/2;
ref.a = ones(32,1)*pi * 0.7;
ref.ph = zeros(32,1);
ncyc = ones(32,2);
tr = 1000;
te = ones(32,1)*10;

t = col(10:10:320);
s = mese(m0, t1, t2, kap, ex, ref, ncyc, tr, te);
figure; 
hold on;
plot(t, col(abs(s)), 'b-');
plot(t, exp(-t./t2)*(1-exp(-tr/t1)), 'r--');
xlabel('time (ms)');
ylabel('signal (a.u.)');
legend('epg','ideal');
hold off;
