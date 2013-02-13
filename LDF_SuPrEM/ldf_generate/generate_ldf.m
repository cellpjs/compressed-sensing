function ldf_mat = generate_ldf(M, N, dv, peg)

% Generates an N x M LDF with regular left degree dv and regular right
% degree dc. Note M > N and M*dv = N*dc
% If it fails to generate an exactly regular matrix in two attempts, 
% outputs an error message.
%
% Parameters:
% - M = Length of sparse vectors to be measured
% - N = Number of observations
% - dv = Variable node degree
% - peg = If this is 1, uses progressive edge growth (PEG) algorithm and avoids
% short loops. If this is 0, randomly generates. PEG algorithm gives
% superior results, however it may not always produce regular LDFs. Random
% method is also faster. Default uses peg = 1.
% Written by Jinsoo Park and Mehmet Akcakaya, 2009.

if (nargin < 4)
    peg = 1;
end

dc = round(M*dv/N);

if (M*dv - N *dc) ~= 0
    error('Invalid parameters: M*dv/N has to be an integer')
end
K = M*dv;

if peg == 1
    disp('Attempting to generate regular LDF using progressive edge growth algorithm...');
    F = peg_for_gen(M,N,dv,dc);

    f1 = sum(F);
    f2 = sum(F');

    if (sum(f1 == f1(1)) ~= N) || (sum(f2 == f2(1)) ~= M)
        disp('Failed to generate regular LDF. Attempting again using progressive edge growth algorithm...')
        disp('Attempt 2 to generate regular LDF...');
        F = peg_for_gen(M,N,dv,dc);
        f1 = sum(F);
        f2 = sum(F');
        if (sum(f1 == f1(1)) ~= N) || (sum(f2 == f2(1)) ~= M)
            error('Failed to generate regular LDF. Please try again later...')
        end
    end
    
elseif peg == 0
    disp('Attempting to generate regular LDF randomly...');
    F = random_for_gen(M,N,dv,dc);

    f1 = sum(F);
    f2 = sum(F');

    if (sum(f1 == f1(1)) ~= N) || (sum(f2 == f2(1)) ~= M)
        disp('Failed to generate regular LDF. Attempting random generation again...')
        disp('Attempt 2 to generate regular LDF...');
        F = random_for_gen(M,N,dv,dc);
        f1 = sum(F);
        f2 = sum(F');
        if (sum(f1 == f1(1)) ~= N) || (sum(f2 == f2(1)) ~= M)
            error('Failed to generate regular LDF. Please try again later...')
        end
    end
else
    error('Invalid parameter for peg, please see help menu.')
end

disp('Success! Creating substructures...');
A = F';

E = AtoE_convert(A, M, N, dv, dc);

ldf_mat = struct('A', A, 'E', E, 'M', M, 'N', N, 'dv', dv, 'dc', dc, 'K', K);