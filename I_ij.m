function u_Iij = I_ij(x,q,i,j, N)
u_Iij = zeros(q*N,1);
u_Iij((i-1)*N+1:i*N) =  x((j-1)*N+1:(j)*N);
end



