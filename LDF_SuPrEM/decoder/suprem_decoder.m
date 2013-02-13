function [xhat xhat_debiased] = suprem_decoder(r, ldf_mat, L, sd, tot_no_iter, no_reweigh, stop_cond, debiasing)

% The decoding function. Calculates an M x 1 sparse estimate xhat, when
% an N x 1 vector, r = Ax + n is measured, and x is sparse. 
% - r: Observed vector.
% - ldf_mat: The matrix structure created by generate_ldf
% - L: The sparsity of x. (i.e. ||x||_0)
% - sd: Standard deviation of the noise n (if n is not Gaussian, this can be
% taken an estimate for ||n||_2/sqrt(N)).
% - tot_no_iter: Total number of iterations. 500 by default, requires less in
% practice.
% - no_reweigh: If reweighing is to be used. Default is 10. 1 gives the
% non-reweighted SuPrEM algorithm. 
% - stop_cond: If 1, a simple stopping condition is used. If 0, no stopping
% condition is used. Default is 0.
% - debiasing: If 1 solves the appropriate least-squares problem at the end to get
% a debiased version. If 0 returns all-zeros vector for xhat_debiased. Default is 1.
% 
% Written by Jinsoo Park and Mehmet Akcakaya, 2009.

if (nargin < 5)
    tot_no_iter = 500;
    no_reweigh = 10;
    stop_cond = 0;
    debiasing = 1;
elseif (nargin < 6)
    no_reweigh = 10;
    stop_cond = 0;
    debiasing = 1;
elseif (nargin < 7)
    stop_cond = 0;
    debiasing = 1;
elseif (nargin < 8)
    debiasing = 1;
end

if length(r) ~= ldf_mat.N
    error('Invalid measurement vector: Please use the same matrix for measurement.');
end

if L ~= round(L)
    error('Sparsity must be an integer')
end

if stop_cond ~= 0 && stop_cond ~= 1
    error('Invalid parameter for stopping condition');
end

if debiasing ~= 0 && debiasing ~= 1
    error('Invalid parameter for debiasing');
end


iterations_per_round = round(tot_no_iter/no_reweigh);
if (iterations_per_round * no_reweigh - tot_no_iter) ~= 0
    error('Invalid parameters: Total number of iterations must be a multiple of the number of reweighings');
end

%disp('Starting Decoder')
time1 = cputime;
xhat = suprem(r, ldf_mat.E, ldf_mat.K, ldf_mat.M, ldf_mat.N, ldf_mat.dc, ldf_mat.dv, L, sd, no_reweigh, iterations_per_round, stop_cond);
time1 = cputime - time1;
%disp(sprintf('Completed in %d seconds', time1));

M = ldf_mat.M;
xhat_debiased = zeros(M, 1);
if debiasing == 1
    F = ldf_mat.A;
    in = (xhat~=0);
    S = F(:,in);
    disp('Debiasing...')
    [xnz flag] = lsqr(S, r, 1e-9, 100);
    xhat_debiased(in) = xnz;
end
