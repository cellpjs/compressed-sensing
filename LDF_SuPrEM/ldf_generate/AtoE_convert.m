function E = AtoE_convert(A, M, N, dv, dc)
      
K = M*dv;
E = zeros(K, 1);
dcount=zeros(1,M);
for i=1:N
    v=find(A(i,:));
    for k=1:dc
        j=v(k);
        vsock=dv*(j-1)+dcount(j);
        dcount(j)=dcount(j)+1;
        E(dc*(i-1)+k)=vsock;
    end
end