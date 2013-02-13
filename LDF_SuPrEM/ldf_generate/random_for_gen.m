function F = random_for_gen(M, N, dl ,dr)

F = speye(M);
F = F(:, 1:N);
for ind0 = 1:N
    F(ind0,ind0) = 0;
end


LH_list = zeros(M, dl);
RH_list = zeros(N, dr + 2);
check_degrees = zeros(1, N);

for ind1 = 1:M
    for k = 1:dl
        
        a = min(check_degrees);
        in1 = find(check_degrees == a);     %find the check nodes with minimum number of edges connected to it
        in_index = ceil(rand(1,1) * length(in1));   %break ties randomly
        check_index = in1(in_index);                %pick the check node

        F(ind1, check_index) = 1;                   %set corresponding element of F to 1
        check_degrees(check_index) = check_degrees(check_index) + 1;    %keep count of number of edges connected to the check node
        LH_list(ind1, k) = check_index;                             %update left hand list
        checks = RH_list(check_index, :);                           %update right hand list by adding ind1
        checks_unique = unique([checks ind1]);
        RH_list(check_index, :) = [checks_unique zeros(1, dr + 2 - length(checks_unique))];
    end
end
