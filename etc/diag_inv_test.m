% script diag_inv_test.m
% test to compare run times of diag_inv versus regular inversion.

% instantiation
max = 1000;
step = 1;
size = [step:step:max];

time_diag = NaN(length(size),1);
time_clas = NaN(length(size),1);

for i = 1:length(size)
    A = randn(size(i));
    S = A'*A;
    
    tic;
    inv_diag = diag_inv(S, 1);
    time_diag(i) = toc*1000;
    
    tic; 
    inv_clas = inv(S);
    inv_clas = inv_clas(1,1);
    time_clas(i) = toc*1000;
    
%     if norm(inv_diag - inv_clas(1,1),2) > 10^-2
%         error('inversions do not match for i=%u!', i);
%     end
end

figure;
plot(size, time_diag, 'b', size, time_clas, 'r');
xlabel('problem dimension');
ylabel('run time (ms)');
legend('schur', 'classical');