% The main test function. Comments are given below.
% Calls two functions: generate_ldf and suprem_decoder
% Written by Jinsoo Park and Mehmet Akcakaya, 2009.


path(path, './ldf_generate');
path(path, './decoder');

%Set LDF parameters
M = 10000; %length of vector
N = 2500;  %number of measurements
dv = 3;    %left degree for LDF

%Generate LDF matrix
ldf_mat = generate_ldf(M, N, dv, 0);
F = ldf_mat.A;

%Set signal and observation parameters
L = 500;    %sparsity level
SNR = 36;
sd = 10^(-SNR/20);  

%Generate signal and (noisy) measurements
x = [randn(L,1); zeros(M-L,1)];
x = x(randperm(M));
r = F*x;
x = x*sqrt(N)/norm(r);  %normalize for SNR
r = F*x;
r = r + sd*randn(N,1);

%Calculate Genie Bound
in =(x~=0);
S = F(:,in);
[xnz flag] = lsqr(S,r,1e-9,100);
xgenie = zeros(M,1);
xgenie(in) = xnz;

%Recovery
%[xhat xhat_debiased] = suprem_decoder(r, ldf_mat, L, sd, 500, 1, 1, 1);  %no reweighing
disp('Starting Decoder')
[xhat xhat_debiased] = suprem_decoder(r, ldf_mat, L, sd); %with reweighing
xdis_norm = norm(xhat-x)/norm(x);
disp(xdis_norm);

if (sum(abs(xhat_debiased) == 0)) == M
    disp('No debiasing was requested. Calculating distortion over genie for xhat...')
    dB_over_genie = 10*log10(norm(xhat-x)^2/norm(xgenie-x)^2);
    display(dB_over_genie);
else
    disp('Calculating distortion over genie for xhat_debiased...')
    dB_over_genie = 10*log10(norm(xhat_debiased-x)^2/norm(xgenie-x)^2);
    display(dB_over_genie);
end

