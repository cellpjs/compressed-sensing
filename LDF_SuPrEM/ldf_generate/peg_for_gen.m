function F = peg_for_gen(M,N,dl,dr);

% clear all;
% M = 200000; 
% dl = 3; dr = 6;
% N = M*dl/dr;

F = speye(M);
F = F(:, 1:N);
for ind0 = 1:N
    F(ind0,ind0) = 0;
end


LH_list = zeros(M, dl);
RH_list = zeros(N, dr + 2);
check_degrees = zeros(1, N);

for ind1 = 1:M
    %if mod(ind1, 20) == 1
    %   ind1
    %end
    for k = 1:dl
        if k == 1;
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
            
            
        else
            n_check = zeros(1, N); %this is the "check node neighborhood" of variable node ind1
            nghd = LH_list(ind1, find(LH_list(ind1, :)));  %these are the neighboring check nodes
            n_check(nghd) = 1;          % add them to the "check node neighborhood"
            n_size = sum(n_check);
            
            logic1 = logical(1);
            logic2 = logical(1);
            logic3 = logical(1);    %for girth 6
            loop_count = 0;
            
            while(logic1 && logic2 && logic3)
                %updates if we ar still going deeper in the graph
                n_check_old = n_check;
                nghd_old = nghd;
                n_size_old = n_size;
                
                %take the var nghd of the previous check nodes
                varnghd_mtx = RH_list(nghd, :);                         %this is a matrix
                varnghd_stack = varnghd_mtx(:);                         %stack them up as a vector
                varnghd = unique(varnghd_stack(varnghd_stack ~= 0));    %take the non-zero elements
                
                %take the check nghd of the var nghd of the previous check nodes
                %we've gone one level deeper in the "check node neighborhood" of variable node ind1
                nghd_mtx = LH_list(varnghd, :);                         
                nghd_stack = nghd_mtx(:);
                nghd = unique(nghd_stack(nghd_stack ~= 0));
                
                n_check(nghd) = 1;
                n_size = sum(n_check);
                
                if n_size < N && n_size == n_size_old
                    logic1 = logical(0);
                end
                
                if n_size_old < N && n_size == N
                    logic2 = logical(0);
                end        
                
                if loop_count > 1
                    logic3 = logical(0);
                end
                loop_count = loop_count + 1;
            end
            
            nghd_comp = 1 - n_check_old;
            original_check_indices = find(nghd_comp);
            checks_in_complement = check_degrees(original_check_indices);
            a = min(checks_in_complement);
            in1 = find(checks_in_complement == a); 
            in_index = ceil(rand(1,1) * length(in1));   %break ties randomly

            check_index = original_check_indices(in1(in_index)); 
            F(ind1, check_index) = 1;                   %set corresponding element of F to 1
            check_degrees(check_index) = check_degrees(check_index) + 1;    %keep count of number of edges connected to the check node
            LH_list(ind1, k) = check_index;                             %update left hand list
            checks = RH_list(check_index, :);                           %update right hand list by adding ind1
            checks_unique = unique([checks ind1]);
            RH_list(check_index, :) = [checks_unique zeros(1, dr + 2 - length(checks_unique))];
         
            
        end
    end
end